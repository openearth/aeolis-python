"""Configuration schema API (dynamic form definition).

Parses ``aeolis/constants.py`` to expose ``DEFAULT_CONFIG`` as a
structured schema: sections (from ``# --- name --- #`` comment lines),
per-parameter defaults, types, units and descriptions (from the inline
``# [unit] description`` comments). Implemented in Phase 2.
"""
