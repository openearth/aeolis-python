/* Colormaps as RGB stop arrays -> 256-entry ramps for WebGL textures
 * and CSS gradients for UI previews. */
"use strict";

const Colormaps = (() => {

  const STOPS = {
    viridis: [[68,1,84],[71,44,122],[59,81,139],[44,113,142],[33,144,141],[39,173,129],[92,200,99],[170,220,50],[253,231,37]],
    turbo: [[48,18,59],[70,107,227],[40,187,235],[74,236,151],[183,246,75],[249,189,38],[245,96,23],[177,25,9],[122,4,3]],
    plasma: [[13,8,135],[84,2,163],[139,10,165],[185,50,137],[219,92,104],[244,136,73],[254,188,43],[240,249,33]],
    RdBu: [[103,0,31],[178,24,43],[214,96,77],[244,165,130],[247,247,247],[146,197,222],[67,147,195],[33,102,172],[5,48,97]],
    gray: [[15,15,15],[240,240,240]],
    terrain: [[40,54,154],[0,120,180],[80,180,120],[210,205,120],[160,110,70],[240,240,240]],
    topo_dutch: [[10,40,90],[40,110,170],[130,190,220],[235,225,180],[210,180,110],[120,160,80],[70,110,50]],
    phase: [[240,80,80],[200,160,40],[80,190,80],[40,170,200],[110,90,220],[220,80,190],[240,80,80]],
    Greens: [[247,252,245],[199,233,192],[161,217,155],[116,196,118],[65,171,93],[35,139,69],[0,90,50]],
    sand: [[255,250,235],[240,220,160],[214,178,110],[181,137,74],[140,98,52],[92,62,34]],
  };

  const cache = new Map();

  function names() { return Object.keys(STOPS); }

  /* A trailing "!r" on a name means "reversed" (used to store the reverse
   * flag inside the cmap string, so every consumer honours it for free). */
  function _resolve(name) {
    let rev = false;
    if (typeof name === "string" && name.endsWith("!r")) { rev = true; name = name.slice(0, -2); }
    const stops = STOPS[name] || STOPS.viridis;
    return rev ? stops.slice().reverse() : stops;
  }

  /* Combine a base colormap name with a reverse flag into a cmap string. */
  function withReverse(name, reversed) {
    const base = String(name || "viridis").replace(/!r$/, "");
    return reversed ? `${base}!r` : base;
  }
  function isReversed(name) { return typeof name === "string" && name.endsWith("!r"); }
  function baseName(name) { return String(name || "viridis").replace(/!r$/, ""); }

  function sample(name, t) {
    const stops = _resolve(name);
    t = U.clamp(t, 0, 1);
    const pos = t * (stops.length - 1);
    const i = Math.min(stops.length - 2, Math.floor(pos));
    const f = pos - i;
    const a = stops[i], b = stops[i + 1];
    return [
      Math.round(a[0] + (b[0] - a[0]) * f),
      Math.round(a[1] + (b[1] - a[1]) * f),
      Math.round(a[2] + (b[2] - a[2]) * f),
    ];
  }

  /* 256x1 RGBA ramp for WebGL upload. */
  function ramp(name) {
    if (cache.has(name)) return cache.get(name);
    const data = new Uint8Array(256 * 4);
    for (let i = 0; i < 256; i += 1) {
      const [r, g, b] = sample(name, i / 255);
      data[i * 4] = r; data[i * 4 + 1] = g; data[i * 4 + 2] = b; data[i * 4 + 3] = 255;
    }
    cache.set(name, data);
    return data;
  }

  function cssGradient(name, direction = "to right") {
    const stops = _resolve(name);
    const parts = stops.map((c, i) =>
      `rgb(${c[0]},${c[1]},${c[2]}) ${(100 * i / (stops.length - 1)).toFixed(1)}%`);
    return `linear-gradient(${direction}, ${parts.join(", ")})`;
  }

  function stops(name) { return _resolve(name); }

  return { names, sample, ramp, cssGradient, stops, withReverse, isReversed, baseName };
})();
