# AeoLiS Web GUI (`aeolis webui`)

The modern AeoLiS GUI: a zero-build web application presented in a
native desktop window. It replaces nothing — the legacy Tkinter GUI
(`aeolis gui`, `aeolis/gui/`) is untouched and still works.

```
pip install -e .[webui]        # adds pywebview (native window)
pip install -e .[webui-data]   # optional: rasterio, pyproj, cdsapi for data downloads
aeolis webui [path/to/aeolis.txt] [--port N] [--browser]
```

Without pywebview the GUI opens in the default browser instead.

## Architecture

**Zero-build frontend.** Vanilla ES2020 JavaScript loaded as plain
`<script>` tags — no Node, no npm, no bundler. The dev loop is: edit a
file under `web/`, press F5. Third-party libraries are vendored as
single files in `web/vendor/` (MapLibre GL for the map, Terra Draw for
polygon drawing, uPlot for graphs, proj4 for client-side CRS
transforms — see `web/vendor/README.md` for versions).

**Stdlib backend.** `backend/httpd.py` runs a threaded
`http.server` with a small route registry; API modules register routes
with the `@route(method, path)` decorator and are imported in
`_load_api_modules()`. JSON for control, raw binary
(`application/octet-stream`) for bulk data (grid meshes, netCDF field
slabs). All netCDF access holds `util.NC_LOCK` (HDF5 is not
thread-safe); the GUI process sets `HDF5_USE_FILE_LOCKING=FALSE` so it
can read output files that a running model is still writing (the run
subprocess gets the variable stripped again).

**Project folder convention.** A project is the folder containing
`aeolis.txt`. Everything the GUI creates lives in one `gui/` subfolder
so a project stays a single portable directory:

```
<project>/
  aeolis.txt  *.grd  wind.txt ...   model inputs
  aeolis.nc  aeolis.log             run output
  gui/
    project.json                    GUI state (CRS, layers, history)
    polygons.json                   objects store (polygons, transects)
    rawdata/  + manifest.json       downloaded raw data + provenance
    cache/                          derived caches (safe to delete)
```

**Single source of truth.** The Settings form is generated at runtime
from `aeolis/constants.py` (`backend/schema_api.py` parses the
`# --- section --- #` headers and `# [unit] description` comments;
defaults come from the evaluated `DEFAULT_CONFIG`), and reading/writing
`aeolis.txt` goes through `aeolis.inout`. New model parameters appear
in the GUI automatically.

**Fast output viewing.** `backend/output_api.py` streams per-timestep
`float32` slabs with an LRU cache; `web/js/core/fieldlayer.js` is a
MapLibre custom WebGL2 layer with ping-pong value buffers and GPU
interpolation (`uFrac`) between timesteps, plus client-side prefetch —
scrubbing the time slider is instantaneous and nothing is pre-rendered,
so zooming stays sharp.

**Coordinate systems.** `web/js/core/crs.js` chains
model CRS ⇄ WGS84 ⇄ WebMercator via proj4. Projected EPSG codes get
basemaps (grey CARTO / Esri satellite); the *local/conceptual* mode
maps model meters directly onto Mercator meters with no basemap, for
academic/fictional cases.

## Layout of the code

```
aeolis/webui/
  launcher.py            port scan, server thread, pywebview window
  backend/
    httpd.py             server + route registry (start here)
    schema_api.py        constants.py -> settings form schema
    config_api.py        aeolis.txt load/save via aeolis.inout
    grd_io.py            grid generate/read/write + shear preview
    grid_api.py          Grid tab endpoints
    domain_api.py        raw data manifest, modify, interpolate
    conditions_api.py    wind/tide/wave generation + downloads
    run_manager.py       RunnerBackend (local subprocess; HPC planned)
    run_api.py           Run tab endpoints
    output_api.py        netCDF/gridfield binary streaming
    project.py/objects…  project folder + shapes store
    datasources/         rws_lidar, jarkus, vaklodingen, xyz_import,
                         tiff_import, synthetic, era5, waterinfo
  web/
    index.html           layout + script load order
    js/core/             state bus, map view, CRS, field layer, draw,
                         panels, playbar, graphs
    js/tabs/             settings, grid, domain, conditions, run, viewer
```

## Adding things

- **A data source**: add a module in `backend/datasources/` with
  `check(bounds, job)` and `download(bounds, years, dest_dir, job)`
  returning manifest entries, and register it in
  `datasources/__init__.py` (`SOURCES` + `INFO`).
- **A run backend** (e.g. the planned Deltares HPC/SLURM runner):
  subclass `RunnerBackend` in `backend/run_manager.py` and add it to
  `BACKENDS`; the Run tab UI picks it up automatically.
- **An API route**: `from aeolis.webui.backend.httpd import route` and
  decorate a handler; import the module in `_load_api_modules()`.

## Tests

```
python -m pytest aeolis/tests/test_webui.py
```
