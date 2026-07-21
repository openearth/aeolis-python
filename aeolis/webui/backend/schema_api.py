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

# hidden from the GUI entirely: orientation is embedded in the grid
# coordinates, nx/ny are derived from the grid files, and output_types
# is covered by the per-variable statistics in the output_vars picker
HIDDEN = {"alfa", "nx", "ny", "output_types"}
READONLY = set()

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

# conditional visibility rules, evaluated live by the frontend.
# A rule is one of:
#   {"key": param, "in": [values]}          value match
#   {"any": [rule, ...]}                    at least one holds
#   {"all": [rule, ...]}                    all hold
_SURF_MOIST = ["fc", "resd_moist", "satw_moist", "satd_moist", "nw_moist",
               "nd_moist", "mw_moist", "md_moist", "alfaw_moist", "alfad_moist",
               "thick_moist"]
_GROUNDWATER = ["boundary_gw", "K_gw", "ne_gw", "D_gw", "tfac_gw", "Cl_gw",
                "in_gw", "GW_stat"]
# parameters of the old (duran) vegetation framework vs the new grass one
_VEG_DURAN = ["avg_time", "gamma_vegshear", "hveg_max", "dzb_opt", "V_ver",
              "germinate", "lateral", "veg_gamma", "okin_c1_veg",
              "okin_initialred_veg", "rhoveg_max", "t_veg", "v_gam"]
_VEG_GRASS = ["veg_res_factor", "dt_veg", "species_names", "d_tiller", "r_stem",
              "alpha_uw", "alpha_Nt", "alpha_0", "G_h", "G_c", "G_s", "Hveg",
              "phi_h", "Nt_max", "R_cov", "lmax_c", "mu_c", "alpha_s", "nu_s",
              "T_burial", "gamma_h", "dzb_tol_c", "dzb_tol_s", "dzb_opt_h",
              "dzb_opt_c", "dzb_opt_s", "beta_veg", "m_veg", "c1_okin",
              "alpha_comp", "T_flood", "gamma_Nt_decay", "pNt_zeta",
              "bounce", "alpha_lift"]

VISIBLE_IF = {
    # vegetation framework: duran uses veg_file, grass uses hveg/Nt
    "veg_file": {"key": "method_vegetation", "in": ["duran"]},
    "hveg_file": {"key": "method_vegetation", "in": ["grass"]},
    "Nt_file": {"key": "method_vegetation", "in": ["grass"]},
    **{k: {"key": "method_vegetation", "in": ["duran"]} for k in _VEG_DURAN},
    **{k: {"key": "method_vegetation", "in": ["grass"]} for k in _VEG_GRASS},
    # sand fence variants of the Okin parameters
    "okin_c1_fence": {"key": "process_fences", "in": [True]},
    "okin_initialred_fence": {"key": "process_fences", "in": [True]},
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
    # boundary flux factors only apply to 'flux' boundaries
    "offshore_flux": {"key": "boundary_offshore", "in": ["flux"]},
    "onshore_flux": {"key": "boundary_onshore", "in": ["flux"]},
    "lateral_flux": {"key": "boundary_lateral", "in": ["flux"]},
    # method selectors follow their process boolean
    "method_shear": {"key": "process_shear", "in": [True]},
    "method_transport": {"key": "process_transport", "in": [True]},
    "method_grainspeed": {"key": "process_transport", "in": [True]},
    "method_moist_process": {"any": [
        {"key": "process_moist", "in": [True]},
        {"key": "th_moisture", "in": [True]}]},
    "method_moist_threshold": {"key": "th_moisture", "in": [True]},
    "method_vegetation": {"key": "process_vegetation", "in": [True]},
    "vegshear_type": {"key": "process_vegetation", "in": [True]},
    "veggrowth_type": {"all": [
        {"key": "process_vegetation", "in": [True]},
        {"key": "method_vegetation", "in": ["duran"]}]},
}

# no whole sections tied to one method anymore (vegetation is merged)
SECTION_VISIBLE_IF = {}

# process switch per section: when the rule is false the section is
# parked below a "Disabled" divider at the bottom of the Settings panel
SECTION_ENABLED_IF = {
    "Topographic steering (shear)": {"key": "process_shear", "in": [True]},
    "Separation bubble": {"key": "process_separation", "in": [True]},
    "Sediment transport": {"key": "process_transport", "in": [True]},
    "Armouring and sheltering": {"key": "th_sheltering", "in": [True]},
    "Hydrodynamics and waves": {"any": [
        {"key": "process_tide", "in": [True]},
        {"key": "process_wave", "in": [True]},
        {"key": "process_runup", "in": [True]},
        {"key": "process_wet_bed_reset", "in": [True]}]},
    "Moisture and groundwater": {"any": [
        {"key": "process_moist", "in": [True]},
        {"key": "th_moisture", "in": [True]},
        {"key": "process_groundwater", "in": [True]}]},
    "Avalanching": {"key": "process_avalanche", "in": [True]},
    "Vegetation": {"key": "process_vegetation", "in": [True]},
    "Bed interaction": {"key": "process_bedinteraction", "in": [True]},
    "Dune erosion": {"key": "process_dune_erosion", "in": [True]},
    "Salt": {"any": [
        {"key": "process_salt", "in": [True]},
        {"key": "th_salt", "in": [True]}]},
}

# per-section docs links removed in favour of one topbar entry point;
# the /api/docs proxy below stays (used by the docs popup)
DOCS_LINKS = {}

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
                "picker": "output_vars" if key == "output_vars" else None,
            })
        if params:
            out_sections.append({
                "name": section["name"],
                "params": params,
                "docs": DOCS_LINKS.get(section["name"]),
                "visible_if": SECTION_VISIBLE_IF.get(section["name"]),
                "enabled_if": SECTION_ENABLED_IF.get(section["name"]),
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
# available output variables (for the output_vars picker)
# ---------------------------------------------------------------------

_VAR_COMMENT_RE = re.compile(r"^\s*'(\w+)',\s*#\s*(.*?)\s*$")
_outvars_cache = None


def list_output_vars():
    """All spatial model-state variables that can be written to the
    netCDF output, with dims and the description comment from
    constants.py. Statistics (avg/sum/var/min/max) can be requested per
    variable with a ``_<stat>`` suffix in output_vars."""
    global _outvars_cache
    if _outvars_cache is not None:
        return _outvars_cache

    descs = {}
    for line in open(aeolis.constants.__file__, "r", encoding="utf-8"):
        match = _VAR_COMMENT_RE.match(line)
        if match and match.group(1) not in descs:
            descs[match.group(1)] = match.group(2)

    out = []
    seen = set()
    for dims, names in aeolis.constants.MODEL_STATE.items():
        if "ny" not in dims or "nx" not in dims:
            continue
        for name in names:
            if name in seen:
                continue
            seen.add(name)
            out.append({
                "name": name,
                "dims": list(dims),
                "desc": descs.get(name, ""),
            })
    out.sort(key=lambda v: v["name"].lower())
    _outvars_cache = out
    return out


@route("GET", "/api/schema/output_vars")
def _output_vars(handler, query, tail):
    send_json(handler, {
        "variables": list_output_vars(),
        "stats": ["avg", "sum", "var", "min", "max"],
    })


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
