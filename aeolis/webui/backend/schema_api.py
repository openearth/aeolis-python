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
}

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
            default = DEFAULT_CONFIG[key]
            info = meta[key]
            params.append({
                "key": key,
                "default": default,
                "type": _infer_type(key, default),
                "unit": info["unit"],
                "desc": info["desc"],
                "options": _infer_options(default, info["desc"]),
                "link": LINKS.get(key),
            })
        if params:
            out_sections.append({"name": section["name"], "params": params})

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
