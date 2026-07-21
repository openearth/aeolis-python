"""Configuration schema API: the dynamic form definition for the
Settings tab.

The schema is derived from ``aeolis/constants.py`` at runtime so the
GUI is always in sync with the code:

- section titles come from the ``# --- name --- #`` comment lines
  inside the ``DEFAULT_CONFIG`` literal,
- descriptions and units come from the inline ``# [unit] description``
  comments,
- default values come from the evaluated ``DEFAULT_CONFIG`` dict itself
  (authoritative, matches Python dict semantics for duplicated keys).

Parameters whose input is produced by another GUI tab (grids, domain
files, boundary-condition timeseries) carry a ``link`` so the form can
render a cross-reference instead of a text input.
"""

import re

import aeolis.constants
from aeolis.constants import DEFAULT_CONFIG
from aeolis.webui.backend.httpd import route
from aeolis.webui.backend.util import send_json

# parameters produced by other tabs
LINKS = {
    "xgrid_file": "grid",
    "ygrid_file": "grid",
    "bed_file": "domain",
    "ne_file": "domain",
    "veg_file": "domain",
    "hveg_file": "domain",
    "Nt_file": "domain",
    "wind_file": "conditions",
    "tide_file": "conditions",
    "wave_file": "conditions",
    "meteo_file": "conditions",
    # 2D spatial files / masks (bedcomp_file is 4D -> stays a text field)
    "threshold_file": "domain",
    "fence_file": "domain",
    "supply_file": "domain",
    "wave_mask": "domain",
    "tide_mask": "domain",
    "runup_mask": "domain",
    "threshold_mask": "domain",
    "gw_mask": "domain",
    "vver_mask": "domain",
}

# hidden from the GUI entirely (orientation is embedded in the grid
# coordinates; nx/ny are derived from the grid files)
HIDDEN = {"alfa"}
READONLY = {"nx", "ny"}

# valid options verified against the model source (see webui README)
OPTIONS_OVERRIDE = {
    # transport.py equilibrium(), lines 382-450
    "method_transport": ["bagnold", "bagnold_gs", "kawamura", "lettau", "dk",
                         "sauermann", "vanrijn_strypsteen"],
    # transport.py grainspeed dispatch
    "method_grainspeed": ["windspeed", "constant", "duran", "duran_full",
                          "duran_uniform", "duran_flat"],
    # hydro.py update(), lines 324-405
    "method_moist_process": ["infiltration", "surf_moisture"],
    # threshold.py compute_moisture(), lines 247-295
    "method_moist_threshold": ["belly_johnson", "hotta", "chepil", "saleh_fryear",
                               "saleh_fryear_mod", "shao", "dong_2002",
                               "gregory_darwish", "cornelis", "dong_2007"],
    # model.py initialize, lines 297-302
    "method_vegetation": ["duran", "grass"],
    # wind.py initialize/shear
    "method_shear": ["fft", "quasi2d", "1Dstacks"],
    # vegetation.py
    "vegshear_type": ["raupach", "okin"],
}

# conditional visibility: param -> {"key": other_param, "in": [values]}
# (evaluated live by the frontend; a section with no visible params is
# hidden as a whole)
_SURF_MOIST = ["fc", "resd_moist", "satw_moist", "satd_moist", "nw_moist",
               "nd_moist", "mw_moist", "md_moist", "alfaw_moist", "alfad_moist",
               "thick_moist"]
_GROUNDWATER = ["boundary_gw", "K_gw", "ne_gw", "D_gw", "tfac_gw", "Cl_gw",
                "in_gw", "GW_stat"]
VISIBLE_IF = {
    # vegetation framework: duran uses veg_file, grass uses hveg/Nt
    "veg_file": {"key": "method_vegetation", "in": ["duran"]},
    "hveg_file": {"key": "method_vegetation", "in": ["grass"]},
    "Nt_file": {"key": "method_vegetation", "in": ["grass"]},
    # transport constants per formulation
    "Cb": {"key": "method_transport", "in": ["bagnold", "bagnold_gs", "vanrijn_strypsteen"]},
    "Ck": {"key": "method_transport", "in": ["kawamura"]},
    "Cl": {"key": "method_transport", "in": ["lettau"]},
    "Cdk": {"key": "method_transport", "in": ["dk"]},
    # shear: only the fft WindShear uses the perturbation params
    "L": {"key": "method_shear", "in": ["fft"]},
    "l": {"key": "method_shear", "in": ["fft"]},
    "buffer_width": {"key": "method_shear", "in": ["fft"]},
    # moisture process methods
    "Tdry": {"key": "method_moist_process", "in": ["infiltration"]},
    **{k: {"key": "method_moist_process", "in": ["surf_moisture"]} for k in _SURF_MOIST},
    **{k: {"key": "process_groundwater", "in": [True]} for k in _GROUNDWATER},
    "w1_5": {"key": "method_moist_threshold",
             "in": ["chepil", "saleh_fryear", "saleh_fryear_mod", "gregory_darwish",
                    "cornelis", "dong_2002", "dong_2007"]},
}

# sections whose whole parameter set belongs to one vegetation method
SECTION_VISIBLE_IF = {
    "Vegetation (OLD)": {"key": "method_vegetation", "in": ["duran"]},
    "Grass vegetation model (new vegetation framework)":
        {"key": "method_vegetation", "in": ["grass"]},
}

# critically selected readthedocs pages per section
_DOCS = "https://aeolis.readthedocs.io/en/update_documentation/user/"
DOCS_LINKS = {
    "Grid files (convention *.grd)": _DOCS + "model_setup.html",
    "Model, grid and time settings": _DOCS + "model_setup.html",
    "Input Timeseries": _DOCS + "model_setup.html",
    "Boundary conditions": _DOCS + "model_setup.html",
    "Output (and coupling) settings": _DOCS + "model_setup.html",
    "Other spatial files / masks": _DOCS + "model_setup.html",
    "Process Booleans (True/False)": _DOCS + "model_description.html",
    "Threshold Booleans (True/False)": _DOCS + "model_description.html",
    "Sediment transport formulations": _DOCS + "model_description.html",
    "Topographic steering (shear)": _DOCS + "model_description.html",
    "Vegetation (OLD)": _DOCS + "model_description.html",
    "Grass vegetation model (new vegetation framework)": _DOCS + "model_description.html",
    "Moisture parameters": _DOCS + "model_description.html",
    "Avalanching": _DOCS + "model_description.html",
    "Hydro and waves": _DOCS + "model_description.html",
}

# time-like parameters that get the date/duration helper tool
TIME_PARAMS = {"tstart", "tstop", "dt", "restart", "output_times", "dzb_interval"}

_UNIT_RE = re.compile(r"^\[([^\]]*)\]\s*")
_PARAM_RE = re.compile(r"^\s*'([^']+)'\s*:")
_SECTION_RE = re.compile(r"^#\s*---+\s*(.*?)\s*-{2,}.*#?\s*$")
# trailing "(a, b or c)" / "(a, b, c)" enumerations in descriptions
_OPTIONS_RE = re.compile(r"\(([^()]+?(?:,| or )[^()]+?)\)\s*$")

_cache = None


def build_schema():
    """Parse constants.py into an ordered section/parameter schema."""
    global _cache
    if _cache is not None:
        return _cache

    lines = open(aeolis.constants.__file__, "r", encoding="utf-8").readlines()

    sections = []          # [{name, params: [key, ...]}]
    meta = {}              # key -> {unit, desc}
    seen = set()
    in_config = False
    current = None

    for line in lines:
        stripped = line.strip()
        if not in_config:
            if "DEFAULT_CONFIG" in line and "=" in line and "{" in line:
                in_config = True
            continue
        if stripped.startswith("}"):
            break

        section_match = _SECTION_RE.match(stripped)
        if stripped.startswith("# ---") and section_match:
            name = section_match.group(1).strip().rstrip("-").strip()
            current = next((s for s in sections if s["name"] == name), None)
            if current is None:
                current = {"name": name, "params": []}
                sections.append(current)
            continue

        param_match = _PARAM_RE.match(line)
        if param_match and not stripped.startswith("#"):
            key = param_match.group(1)
            if key not in DEFAULT_CONFIG:
                continue  # continuation-line token, not a parameter
            comment = ""
            if "#" in line:
                comment = line.split("#", 1)[1].strip()
            if key in seen:
                # duplicated key in constants.py: keep first position,
                # refresh the comment if the earlier one was empty
                if comment and not meta[key]["raw"]:
                    meta[key] = _parse_comment(comment)
                continue
            seen.add(key)
            meta[key] = _parse_comment(comment)
            if current is None:
                current = {"name": "General", "params": []}
                sections.append(current)
            current["params"].append(key)

    # any DEFAULT_CONFIG keys not found in the text (safety net)
    missing = [k for k in DEFAULT_CONFIG if k not in seen]
    if missing:
        sections.append({"name": "Other", "params": missing})
        for key in missing:
            meta[key] = _parse_comment("")

    out_sections = []
    for section in sections:
        params = []
        for key in section["params"]:
            if key in HIDDEN:
                continue
            default = DEFAULT_CONFIG[key]
            info = meta[key]
            params.append({
                "key": key,
                "default": default,
                "type": _infer_type(key, default),
                "unit": info["unit"],
                "desc": info["desc"],
                "options": OPTIONS_OVERRIDE.get(key) or _infer_options(default, info["desc"]),
                "link": LINKS.get(key),
                "readonly": key in READONLY,
                "visible_if": VISIBLE_IF.get(key),
                "time_tool": key in TIME_PARAMS,
            })
        if params:
            out_sections.append({
                "name": section["name"],
                "params": params,
                "docs": DOCS_LINKS.get(section["name"]),
                "visible_if": SECTION_VISIBLE_IF.get(section["name"]),
            })

    _cache = {"sections": out_sections}
    return _cache


def _parse_comment(comment):
    raw = comment
    unit = None
    match = _UNIT_RE.match(comment)
    if match:
        unit = match.group(1).strip() or None
        comment = comment[match.end():]
    return {"unit": unit, "desc": comment.strip(), "raw": raw}


def _infer_type(key, default):
    if isinstance(default, bool):
        return "bool"
    if isinstance(default, int):
        return "int"
    if isinstance(default, float):
        return "float"
    if isinstance(default, (list, tuple)):
        return "list"
    if key.endswith(("_file", "_mask")):
        return "file"
    if isinstance(default, str):
        return "str"
    return "any"   # None without a file hint


def _infer_options(default, desc):
    """Extract '(a, b or c)' enumerations for string parameters."""
    if not isinstance(default, str):
        return None
    match = _OPTIONS_RE.search(desc)
    if not match:
        return None
    parts = re.split(r",| or ", match.group(1))
    options = [p.strip() for p in parts if p.strip()]
    if len(options) < 2 or any(" " in o or len(o) > 24 for o in options):
        return None
    if default not in options:
        options.insert(0, default)
    return options


@route("GET", "/api/schema")
def _schema(handler, query, tail):
    send_json(handler, build_schema())


# ---------------------------------------------------------------------
# Documentation proxy: serve readthedocs pages same-origin so the
# floating docs popup can iframe them (RTD error pages set
# X-Frame-Options), with an in-memory cache to avoid rate limits.
# ---------------------------------------------------------------------

_DOCS_CACHE = {}
_ALLOWED_DOCS_HOST = "aeolis.readthedocs.io"


@route("GET", "/api/docs")
def _docs_proxy(handler, query, tail):
    import urllib.parse
    import urllib.request

    from aeolis.webui.backend.util import send_bytes, send_error_json

    url = query.get("url", "")
    parsed = urllib.parse.urlparse(url)
    if parsed.scheme != "https" or parsed.hostname != _ALLOWED_DOCS_HOST:
        send_error_json(handler, "only aeolis.readthedocs.io pages can be shown", 403)
        return

    html = _DOCS_CACHE.get(url)
    if html is None:
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "aeolis-webui"})
            with urllib.request.urlopen(req, timeout=20) as response:
                html = response.read()
        except Exception as exc:  # noqa: BLE001 - shown in the popup
            send_error_json(handler, f"could not load documentation: {exc}", 502)
            return
        # make relative assets/links resolve against the original site
        base = f'<base href="{url}">'.encode()
        lowered = html.lower()
        head = lowered.find(b"<head")
        if head != -1:
            insert = lowered.find(b">", head) + 1
            html = html[:insert] + base + html[insert:]
        _DOCS_CACHE[url] = html
        while len(_DOCS_CACHE) > 30:
            _DOCS_CACHE.pop(next(iter(_DOCS_CACHE)))

    send_bytes(handler, html, ctype="text/html; charset=utf-8")
