/* FieldLayer: MapLibre custom WebGL2 layer rendering a (curvilinear)
 * grid of scalar values with a colormap — used for model output,
 * domain .grd files and raw rasters.
 *
 * Fast scrubbing recipe (ported from the SedTRAILS GUI): two value
 * buffers (ping-pong) with GPU interpolation via uFrac, so animating
 * between loaded timesteps is free; per-timestep slabs stream as raw
 * float32 and are cached/prefetched by the caller. No pre-rendering:
 * the field is drawn from the actual data at every zoom level.
 *
 * NaN cells carry the sentinel -1e30 and are discarded in the shader.
 */
"use strict";

const FieldLayer = (() => {

  const SENTINEL_CUT = -1.0e29;

  const VS = `#version 300 es
    precision highp float;
    uniform mat4 uMatrix;
    uniform float uFrac;
    in vec2 aPos;
    in float aValA;
    in float aValB;
    out float vVal;
    void main() {
      float a = aValA;
      float b = aValB;
      // if either frame is missing data, use the other
      float v = (a < ${SENTINEL_CUT}) ? b
              : (b < ${SENTINEL_CUT}) ? a
              : mix(a, b, uFrac);
      vVal = v;
      gl_Position = uMatrix * vec4(aPos, 0.0, 1.0);
    }`;

  const FS = `#version 300 es
    precision highp float;
    uniform sampler2D uRamp;
    uniform float uMin;
    uniform float uMax;
    uniform float uOpacity;
    in float vVal;
    out vec4 outColor;
    void main() {
      if (vVal < ${SENTINEL_CUT}) discard;
      float t = clamp((vVal - uMin) / max(uMax - uMin, 1e-12), 0.0, 1.0);
      vec4 c = texture(uRamp, vec2(t, 0.5));
      outColor = vec4(c.rgb * uOpacity, uOpacity);
    }`;

  class Layer {
    /* mesh: {x: Float32Array, y: Float32Array, n, s} in model CRS */
    constructor(id, mesh, style = {}) {
      this.id = id;
      this.type = "custom";
      this.renderingMode = "2d";
      this.mesh = mesh;
      this.style = { cmap: "viridis", min: 0, max: 1, opacity: 1, ...style };
      this.visible = true;
      this._pending = { a: null, b: null, frac: 0 };
      this._gl = null;
    }

    /* ---- MapLibre hooks ---- */

    onAdd(map, gl) {
      this._gl = gl;
      this.map = map;
      const compile = (type, src) => {
        const sh = gl.createShader(type);
        gl.shaderSource(sh, src);
        gl.compileShader(sh);
        if (!gl.getShaderParameter(sh, gl.COMPILE_STATUS)) {
          throw new Error(gl.getShaderInfoLog(sh));
        }
        return sh;
      };
      this.prog = gl.createProgram();
      gl.attachShader(this.prog, compile(gl.VERTEX_SHADER, VS));
      gl.attachShader(this.prog, compile(gl.FRAGMENT_SHADER, FS));
      gl.linkProgram(this.prog);
      if (!gl.getProgramParameter(this.prog, gl.LINK_STATUS)) {
        throw new Error(gl.getProgramInfoLog(this.prog));
      }
      this.loc = {
        matrix: gl.getUniformLocation(this.prog, "uMatrix"),
        frac: gl.getUniformLocation(this.prog, "uFrac"),
        ramp: gl.getUniformLocation(this.prog, "uRamp"),
        min: gl.getUniformLocation(this.prog, "uMin"),
        max: gl.getUniformLocation(this.prog, "uMax"),
        opacity: gl.getUniformLocation(this.prog, "uOpacity"),
        pos: gl.getAttribLocation(this.prog, "aPos"),
        valA: gl.getAttribLocation(this.prog, "aValA"),
        valB: gl.getAttribLocation(this.prog, "aValB"),
      };

      this._buildGeometry(gl);

      this.bufValA = gl.createBuffer();
      this.bufValB = gl.createBuffer();
      const nVerts = this.mesh.n * this.mesh.s;
      const zero = new Float32Array(nVerts).fill(-1e30);
      for (const buf of [this.bufValA, this.bufValB]) {
        gl.bindBuffer(gl.ARRAY_BUFFER, buf);
        gl.bufferData(gl.ARRAY_BUFFER, zero, gl.DYNAMIC_DRAW);
      }

      this.rampTex = gl.createTexture();
      this._uploadRamp(gl, this.style.cmap);

      if (this._pending.a) this.setFrames(this._pending.a, this._pending.b, this._pending.frac);
    }

    _buildGeometry(gl) {
      const { x, y, n, s } = this.mesh;
      const pos = new Float32Array(n * s * 2);
      for (let i = 0; i < n * s; i += 1) {
        const ll = CRS.toLngLat([x[i], y[i]]);
        const mc = maplibregl.MercatorCoordinate.fromLngLat({ lng: ll[0], lat: ll[1] });
        pos[i * 2] = mc.x;
        pos[i * 2 + 1] = mc.y;
      }
      this.bufPos = gl.createBuffer();
      gl.bindBuffer(gl.ARRAY_BUFFER, this.bufPos);
      gl.bufferData(gl.ARRAY_BUFFER, pos, gl.STATIC_DRAW);

      // triangle indices (2 per cell)
      const idx = new Uint32Array((n - 1) * (s - 1) * 6);
      let k = 0;
      for (let j = 0; j < n - 1; j += 1) {
        for (let i = 0; i < s - 1; i += 1) {
          const v0 = j * s + i;
          idx[k++] = v0; idx[k++] = v0 + 1; idx[k++] = v0 + s;
          idx[k++] = v0 + 1; idx[k++] = v0 + s + 1; idx[k++] = v0 + s;
        }
      }
      this.nIndices = idx.length;
      this.bufIdx = gl.createBuffer();
      gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, this.bufIdx);
      gl.bufferData(gl.ELEMENT_ARRAY_BUFFER, idx, gl.STATIC_DRAW);
    }

    _uploadRamp(gl, cmap) {
      gl.bindTexture(gl.TEXTURE_2D, this.rampTex);
      gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA, 256, 1, 0, gl.RGBA,
        gl.UNSIGNED_BYTE, Colormaps.ramp(cmap));
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.LINEAR);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.LINEAR);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
    }

    onRemove(map, gl) {
      for (const buf of [this.bufPos, this.bufIdx, this.bufValA, this.bufValB]) {
        if (buf) gl.deleteBuffer(buf);
      }
      if (this.rampTex) gl.deleteTexture(this.rampTex);
      if (this.prog) gl.deleteProgram(this.prog);
      this._gl = null;
    }

    render(gl, args) {
      if (!this.visible || !this.prog) return;
      const matrix = args && args.defaultProjectionData
        ? args.defaultProjectionData.mainMatrix
        : args;

      gl.useProgram(this.prog);
      gl.uniformMatrix4fv(this.loc.matrix, false, matrix);
      gl.uniform1f(this.loc.frac, this._pending.frac || 0);
      gl.uniform1f(this.loc.min, this.style.min);
      gl.uniform1f(this.loc.max, this.style.max);
      gl.uniform1f(this.loc.opacity, this.style.opacity);

      gl.activeTexture(gl.TEXTURE0);
      gl.bindTexture(gl.TEXTURE_2D, this.rampTex);
      gl.uniform1i(this.loc.ramp, 0);

      gl.bindBuffer(gl.ARRAY_BUFFER, this.bufPos);
      gl.enableVertexAttribArray(this.loc.pos);
      gl.vertexAttribPointer(this.loc.pos, 2, gl.FLOAT, false, 0, 0);

      gl.bindBuffer(gl.ARRAY_BUFFER, this.bufValA);
      gl.enableVertexAttribArray(this.loc.valA);
      gl.vertexAttribPointer(this.loc.valA, 1, gl.FLOAT, false, 0, 0);

      gl.bindBuffer(gl.ARRAY_BUFFER, this.bufValB);
      gl.enableVertexAttribArray(this.loc.valB);
      gl.vertexAttribPointer(this.loc.valB, 1, gl.FLOAT, false, 0, 0);

      gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, this.bufIdx);
      gl.enable(gl.BLEND);
      gl.blendFunc(gl.ONE, gl.ONE_MINUS_SRC_ALPHA);
      gl.drawElements(gl.TRIANGLES, this.nIndices, gl.UNSIGNED_INT, 0);
    }

    /* ---- data & style updates ---- */

    /* frameA/frameB: Float32Array of n*s values; frac in [0,1]. */
    setFrames(frameA, frameB, frac = 0) {
      this._pending = { a: frameA, b: frameB || frameA, frac };
      const gl = this._gl;
      if (!gl) return;
      gl.bindBuffer(gl.ARRAY_BUFFER, this.bufValA);
      gl.bufferData(gl.ARRAY_BUFFER, frameA, gl.DYNAMIC_DRAW);
      gl.bindBuffer(gl.ARRAY_BUFFER, this.bufValB);
      gl.bufferData(gl.ARRAY_BUFFER, frameB || frameA, gl.DYNAMIC_DRAW);
      this.repaint();
    }

    setFrac(frac) {
      this._pending.frac = frac;
      this.repaint();
    }

    setStyle(patch) {
      Object.assign(this.style, patch);
      if (patch.cmap && this._gl) this._uploadRamp(this._gl, patch.cmap);
      this.repaint();
    }

    setVisible(visible) {
      this.visible = visible;
      this.repaint();
    }

    repaint() {
      if (this.map) this.map.triggerRepaint();
    }
  }

  /* ---- registry on the map ---- */

  const layers = new Map();

  function create(id, mesh, style, beforeId = null) {
    remove(id);
    const layer = new Layer(id, mesh, style);
    const map = MapView.instance();
    const before = beforeId && map.getLayer(beforeId) ? beforeId : undefined;
    map.addLayer(layer, before);
    layers.set(id, layer);
    return layer;
  }

  function get(id) { return layers.get(id) || null; }

  function remove(id) {
    const layer = layers.get(id);
    if (layer) {
      const map = MapView.instance();
      if (map.getLayer(id)) map.removeLayer(id);
      layers.delete(id);
    }
  }

  /* Parse a binary gridfield payload -> {mesh, values, range}. */
  function parseGridfield(buffer, headers) {
    const head = new Uint32Array(buffer, 0, 2);
    const n = head[0], s = head[1];
    const count = n * s;
    const x = new Float32Array(buffer, 8, count);
    const y = new Float32Array(buffer, 8 + count * 4, count);
    const v = new Float32Array(buffer, 8 + count * 8, count);
    let range = [0, 1];
    const rangeHeader = headers && headers.get("X-Data-Range");
    if (rangeHeader) range = rangeHeader.split(",").map(Number);
    return { mesh: { x, y, n, s }, values: v, range };
  }

  return { create, get, remove, parseGridfield };
})();
