# Vendored frontend libraries

The AeoLiS web GUI is deliberately zero-build: these libraries are
committed verbatim (single-file distributions) and loaded as plain
`<script>` tags. To upgrade, download the new pinned files from unpkg,
update this table, and test.

| File | Package | Version | License |
|------|---------|---------|---------|
| `maplibre-gl.js` / `maplibre-gl.css` | [maplibre-gl](https://github.com/maplibre/maplibre-gl-js) | 5.24.0 | BSD-3-Clause |
| `terra-draw.umd.js` | [terra-draw](https://github.com/JamesLMilner/terra-draw) | 1.32.1 | MIT |
| `terra-draw-maplibre-gl-adapter.umd.js` | [terra-draw-maplibre-gl-adapter](https://github.com/JamesLMilner/terra-draw) | 1.4.1 | MIT |
| `uplot.iife.min.js` / `uplot.min.css` | [uPlot](https://github.com/leeoniya/uPlot) | 1.6.32 | MIT |
| `proj4.js` | [proj4js](https://github.com/proj4js/proj4js) | 2.20.9 | MIT |

Load order matters: `terra-draw.umd.js` must load before the MapLibre
adapter; `maplibre-gl.js` before both.
