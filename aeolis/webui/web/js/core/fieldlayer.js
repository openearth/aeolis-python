/* FieldLayer: MapLibre custom WebGL2 layer rendering a (curvilinear)
 * grid of scalar values with a colormap — used for model output,
 * domain .grd files and raw rasters.
 *
 * Rendering is pcolormesh-style ("nearest" shading): every data point
 * owns a flat-coloured quad — its model cell, with edges halfway to the
 * neighbouring points and a half-cell rim around the border — so the
 * map shows the true model structure instead of values interpolated
 * between points. NaN cells become crisp rectangular holes.
 *
 * Values live in two R32F textures (ping-pong) mixed by uFrac in the
 * shader, so animating between loaded timesteps is free; per-timestep
 * slabs stream as raw float32 and are cached/prefetched by the caller.
 * No pre-rendering: the field is drawn from the actual data at every
 * zoom level.
 *
 * style.mode = "dots" draws round GPU point sprites at the data points
 * instead of cells (scatter view for sample data); style.dotSize sets
 * their diameter in CSS pixels.
 *
 * NaN cells carry the sentinel -1e30 and are discarded in the shader.
 */
"use strict";

const FieldLayer = (() => {

  const SENTINEL_CUT = -1.0e29;

  /* value lookup shared by both fragment paths */
  const FS_COLOR = `
    uniform sampler2D uRamp;
    uniform float uMin;
    uniform float uMax;
    uniform float uOpacity;
    vec4 rampColor(float v) {
      float t = clamp((v - uMin) / max(uMax - uMin, 1e-12), 0.0, 1.0);
      vec4 c = texture(uRamp, vec2(t, 0.5));
      return vec4(c.rgb * uOpacity, uOpacity);
    }`;

  const VS_CELLS = `#version 300 es
    precision highp float;
    uniform mat4 uMatrix;
    uniform int uCols;             // corner columns = s + 1
    in vec2 aPos;
    out vec2 vUv;                  // corner index (i, j) — cell = floor(uv)
    void main() {
      int i = gl_VertexID % uCols;
      int j = gl_VertexID / uCols;
      vUv = vec2(float(i), float(j));
      gl_Position = uMatrix * vec4(aPos, 0.0, 1.0);
    }`;

  const FS_CELLS = `#version 300 es
    precision highp float;
    uniform sampler2D uTexA;
    uniform sampler2D uTexB;
    uniform float uFrac;
    in vec2 vUv;
    out vec4 outColor;
    ${FS_COLOR}
    void main() {
      ivec2 sz = textureSize(uTexA, 0);
      ivec2 cell = clamp(ivec2(floor(vUv)), ivec2(0), sz - 1);
      float a = texelFetch(uTexA, cell, 0).r;
      float b = texelFetch(uTexB, cell, 0).r;
      // if either frame is missing data, use the other
      float v = (a < ${SENTINEL_CUT}) ? b
              : (b < ${SENTINEL_CUT}) ? a
              : mix(a, b, uFrac);
      if (v < ${SENTINEL_CUT}) discard;
      outColor = rampColor(v);
    }`;

  const VS_DOTS = `#version 300 es
    precision highp float;
    uniform mat4 uMatrix;
    uniform sampler2D uTexA;
    uniform sampler2D uTexB;
    uniform float uFrac;
    uniform float uSize;
    in vec2 aPos;
    flat out float vVal;
    void main() {
      ivec2 sz = textureSize(uTexA, 0);
      ivec2 cell = ivec2(gl_VertexID % sz.x, gl_VertexID / sz.x);
      float a = texelFetch(uTexA, cell, 0).r;
      float b = texelFetch(uTexB, cell, 0).r;
      vVal = (a < ${SENTINEL_CUT}) ? b
           : (b < ${SENTINEL_CUT}) ? a
           : mix(a, b, uFrac);
      gl_Position = uMatrix * vec4(aPos, 0.0, 1.0);
      gl_PointSize = uSize;
    }`;

  const FS_DOTS = `#version 300 es
    precision highp float;
    flat in float vVal;
    out vec4 outColor;
    ${FS_COLOR}
    void main() {
      if (vVal < ${SENTINEL_CUT}) discard;
      vec2 d = gl_PointCoord - 0.5;      // round sprite
      if (dot(d, d) > 0.25) discard;
      outColor = rampColor(vVal);
    }`;

  class Layer {
    /* mesh: {x: Float32Array, y: Float32Array, n, s} in model CRS */
    constructor(id, mesh, style = {}) {
      this.id = id;
      this.type = "custom";
      this.renderingMode = "2d";
      this.mesh = mesh;
      this.style = { cmap: "viridis", min: 0, max: 1, opacity: 1,
        mode: "cells", dotSize: 6, ...style };
      this.visible = true;
      this._pending = { a: null, b: null, frac: 0 };
      this._gl = null;
    }

    /* ---- MapLibre hooks ---- */

    onAdd(map, gl) {
      this._gl = gl;
      this.map = map;
      this._cells = this._makeProgram(gl, VS_CELLS, FS_CELLS,
        ["uMatrix", "uCols", "uTexA", "uTexB", "uFrac", "uRamp", "uMin", "uMax", "uOpacity"]);
      this._dots = this._makeProgram(gl, VS_DOTS, FS_DOTS,
        ["uMatrix", "uTexA", "uTexB", "uFrac", "uSize", "uRamp", "uMin", "uMax", "uOpacity"]);

      this._buildCellGeometry(gl);

      // ping-pong value textures (R32F, one texel per data point)
      const { n, s } = this.mesh;
      this.texA = this._makeValueTex(gl, s, n);
      this.texB = this._makeValueTex(gl, s, n);

      this.rampTex = gl.createTexture();
      this._uploadRamp(gl, this.style.cmap);

      if (this._pending.a) this.setFrames(this._pending.a, this._pending.b, this._pending.frac);
    }

    _makeProgram(gl, vsSrc, fsSrc, uniforms) {
      const compile = (type, src) => {
        const sh = gl.createShader(type);
        gl.shaderSource(sh, src);
        gl.compileShader(sh);
        if (!gl.getShaderParameter(sh, gl.COMPILE_STATUS)) {
          throw new Error(gl.getShaderInfoLog(sh));
        }
        return sh;
      };
      const prog = gl.createProgram();
      gl.attachShader(prog, compile(gl.VERTEX_SHADER, vsSrc));
      gl.attachShader(prog, compile(gl.FRAGMENT_SHADER, fsSrc));
      gl.linkProgram(prog);
      if (!gl.getProgramParameter(prog, gl.LINK_STATUS)) {
        throw new Error(gl.getProgramInfoLog(prog));
      }
      const loc = { pos: gl.getAttribLocation(prog, "aPos") };
      for (const u of uniforms) loc[u] = gl.getUniformLocation(prog, u);
      return { prog, loc };
    }

    /* Mercator anchor at the mesh centre. Vertex positions are stored
     * RELATIVE to this point: absolute web-mercator coordinates (~0.5)
     * quantize to ~1 m in a float32 buffer, which made small-dx cells
     * wiggle and crawl while zooming. Relative coords keep sub-mm
     * precision; the anchor goes back in via the (float64) matrix. */
    _anchorMerc() {
      if (!this._anchor) {
        const { x, y, n, s } = this.mesh;
        const mid = (n >> 1) * s + (s >> 1);
        const ll = CRS.toLngLat([x[mid], y[mid]]);
        const mc = maplibregl.MercatorCoordinate.fromLngLat({ lng: ll[0], lat: ll[1] });
        this._anchor = [mc.x, mc.y];
      }
      return this._anchor;
    }

    /* uMatrix' = uMatrix * translate(anchor), folded in float64 JS so
     * the GPU only ever sees small relative vertex coordinates. */
    _shiftedMatrix(matrix) {
      const [ax, ay] = this._anchorMerc();
      const m = Array.from(matrix);
      m[12] = matrix[0] * ax + matrix[4] * ay + matrix[12];
      m[13] = matrix[1] * ax + matrix[5] * ay + matrix[13];
      m[14] = matrix[2] * ax + matrix[6] * ay + matrix[14];
      m[15] = matrix[3] * ax + matrix[7] * ay + matrix[15];
      return m;
    }

    /* Extrapolated centre coordinate at (j, i) with j in [-1..n], i in
     * [-1..s]: linear extension of the outermost spacing, so the border
     * cells get their half-cell rim just like pcolormesh. */
    _ext(arr, j, i) {
      const { n, s } = this.mesh;
      if (j < 0) return 2 * this._ext(arr, 0, i) - this._ext(arr, Math.min(1, n - 1), i);
      if (j >= n) return 2 * this._ext(arr, n - 1, i) - this._ext(arr, Math.max(n - 2, 0), i);
      if (i < 0) return 2 * arr[j * s] - arr[j * s + Math.min(1, s - 1)];
      if (i >= s) return 2 * arr[j * s + s - 1] - arr[j * s + Math.max(s - 2, 0)];
      return arr[j * s + i];
    }

    /* Corner grid: (n+1) x (s+1) cell-corner vertices, each the average
     * of its 4 surrounding (extrapolated) data points; 2 triangles per
     * data cell. The flat cell value is fetched per fragment from the
     * value texture via floor(corner-index uv). */
    _buildCellGeometry(gl) {
      const { x, y, n, s } = this.mesh;
      const cols = s + 1, rows = n + 1;
      const pos = new Float32Array(rows * cols * 2);
      const [ax, ay] = this._anchorMerc();
      let k = 0;
      for (let j = 0; j < rows; j += 1) {
        for (let i = 0; i < cols; i += 1) {
          const cx = (this._ext(x, j - 1, i - 1) + this._ext(x, j - 1, i)
                    + this._ext(x, j, i - 1) + this._ext(x, j, i)) / 4;
          const cy = (this._ext(y, j - 1, i - 1) + this._ext(y, j - 1, i)
                    + this._ext(y, j, i - 1) + this._ext(y, j, i)) / 4;
          const ll = CRS.toLngLat([cx, cy]);
          const mc = maplibregl.MercatorCoordinate.fromLngLat({ lng: ll[0], lat: ll[1] });
          pos[k++] = mc.x - ax;
          pos[k++] = mc.y - ay;
        }
      }
      this.bufCorners = gl.createBuffer();
      gl.bindBuffer(gl.ARRAY_BUFFER, this.bufCorners);
      gl.bufferData(gl.ARRAY_BUFFER, pos, gl.STATIC_DRAW);

      const idx = new Uint32Array(n * s * 6);
      k = 0;
      for (let j = 0; j < n; j += 1) {
        for (let i = 0; i < s; i += 1) {
          const v0 = j * cols + i;
          idx[k++] = v0; idx[k++] = v0 + 1; idx[k++] = v0 + cols;
          idx[k++] = v0 + 1; idx[k++] = v0 + cols + 1; idx[k++] = v0 + cols;
        }
      }
      this.nIndices = idx.length;
      this.bufIdx = gl.createBuffer();
      gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, this.bufIdx);
      gl.bufferData(gl.ELEMENT_ARRAY_BUFFER, idx, gl.STATIC_DRAW);
    }

    /* Data-point positions for dots mode (built on first use). */
    _ensureDotGeometry(gl) {
      if (this.bufCenters) return;
      const { x, y, n, s } = this.mesh;
      const pos = new Float32Array(n * s * 2);
      const [ax, ay] = this._anchorMerc();
      for (let i = 0; i < n * s; i += 1) {
        const ll = CRS.toLngLat([x[i], y[i]]);
        const mc = maplibregl.MercatorCoordinate.fromLngLat({ lng: ll[0], lat: ll[1] });
        pos[i * 2] = mc.x - ax;
        pos[i * 2 + 1] = mc.y - ay;
      }
      this.bufCenters = gl.createBuffer();
      gl.bindBuffer(gl.ARRAY_BUFFER, this.bufCenters);
      gl.bufferData(gl.ARRAY_BUFFER, pos, gl.STATIC_DRAW);
    }

    _makeValueTex(gl, w, h) {
      const tex = gl.createTexture();
      gl.bindTexture(gl.TEXTURE_2D, tex);
      const fill = new Float32Array(w * h).fill(-1e30);
      gl.texImage2D(gl.TEXTURE_2D, 0, gl.R32F, w, h, 0, gl.RED, gl.FLOAT, fill);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
      return tex;
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
      for (const buf of [this.bufCorners, this.bufCenters, this.bufIdx]) {
        if (buf) gl.deleteBuffer(buf);
      }
      for (const tex of [this.texA, this.texB, this.rampTex]) {
        if (tex) gl.deleteTexture(tex);
      }
      for (const p of [this._cells, this._dots]) {
        if (p && p.prog) gl.deleteProgram(p.prog);
      }
      this.bufCenters = null;
      this._gl = null;
    }

    _bindCommon(gl, loc) {
      gl.uniform1f(loc.uFrac, this._pending.frac || 0);
      gl.uniform1f(loc.uMin, this.style.min);
      gl.uniform1f(loc.uMax, this.style.max);
      gl.uniform1f(loc.uOpacity, this.style.opacity);
      gl.activeTexture(gl.TEXTURE0);
      gl.bindTexture(gl.TEXTURE_2D, this.rampTex);
      gl.uniform1i(loc.uRamp, 0);
      gl.activeTexture(gl.TEXTURE1);
      gl.bindTexture(gl.TEXTURE_2D, this.texA);
      gl.uniform1i(loc.uTexA, 1);
      gl.activeTexture(gl.TEXTURE2);
      gl.bindTexture(gl.TEXTURE_2D, this.texB);
      gl.uniform1i(loc.uTexB, 2);
      gl.enable(gl.BLEND);
      gl.blendFunc(gl.ONE, gl.ONE_MINUS_SRC_ALPHA);
    }

    render(gl, args) {
      if (!this.visible || !this._cells) return;
      const matrix = this._shiftedMatrix(args && args.defaultProjectionData
        ? args.defaultProjectionData.mainMatrix
        : args);

      if (this.style.mode === "dots") {
        this._ensureDotGeometry(gl);
        const { prog, loc } = this._dots;
        gl.useProgram(prog);
        gl.uniformMatrix4fv(loc.uMatrix, false, matrix);
        gl.uniform1f(loc.uSize, this.style.dotSize * (window.devicePixelRatio || 1));
        this._bindCommon(gl, loc);
        gl.bindBuffer(gl.ARRAY_BUFFER, this.bufCenters);
        gl.enableVertexAttribArray(loc.pos);
        gl.vertexAttribPointer(loc.pos, 2, gl.FLOAT, false, 0, 0);
        gl.drawArrays(gl.POINTS, 0, this.mesh.n * this.mesh.s);
        return;
      }

      const { prog, loc } = this._cells;
      gl.useProgram(prog);
      gl.uniformMatrix4fv(loc.uMatrix, false, matrix);
      gl.uniform1i(loc.uCols, this.mesh.s + 1);
      this._bindCommon(gl, loc);
      gl.bindBuffer(gl.ARRAY_BUFFER, this.bufCorners);
      gl.enableVertexAttribArray(loc.pos);
      gl.vertexAttribPointer(loc.pos, 2, gl.FLOAT, false, 0, 0);
      gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, this.bufIdx);
      gl.drawElements(gl.TRIANGLES, this.nIndices, gl.UNSIGNED_INT, 0);
    }

    /* ---- data & style updates ---- */

    /* frameA/frameB: Float32Array of n*s values; frac in [0,1]. */
    setFrames(frameA, frameB, frac = 0) {
      this._pending = { a: frameA, b: frameB || frameA, frac };
      const gl = this._gl;
      if (!gl) return;
      const { n, s } = this.mesh;
      gl.bindTexture(gl.TEXTURE_2D, this.texA);
      gl.texSubImage2D(gl.TEXTURE_2D, 0, 0, 0, s, n, gl.RED, gl.FLOAT, frameA);
      gl.bindTexture(gl.TEXTURE_2D, this.texB);
      gl.texSubImage2D(gl.TEXTURE_2D, 0, 0, 0, s, n, gl.RED, gl.FLOAT, frameB || frameA);
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
