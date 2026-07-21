"""Synthetic boundary-condition generators.

Produce AeoLiS-compatible time series arrays for wind, water levels and
waves. All generators share the same time base: seconds relative to the
configuration reference date, spanning [tstart, tstop] at interval dt.

Profile types
-------------
scalar signals (wind speed, water level, Hs, Tp):
    constant  {value}
    blocks    {values: [...], block_duration}
    harmonic  {mean, amplitude, period, phase?}
    linear    {start, end}
wind direction adds:
    rotational {start, rate}      # deg, deg per hour (full rotations)
"""

import numpy as np


def time_axis(tstart, tstop, dt):
    if tstop <= tstart:
        raise ValueError("tstop must be greater than tstart")
    if dt <= 0:
        raise ValueError("dt must be positive")
    n = int(np.floor((tstop - tstart) / dt)) + 1
    if n > 2_000_000:
        raise ValueError(f"time series too long ({n} steps) - increase dt")
    return tstart + np.arange(n) * dt


def profile(t, spec):
    """Evaluate a profile spec on time axis *t* (seconds)."""
    kind = spec.get("type", "constant")
    p = spec

    if kind == "constant":
        return np.full_like(t, float(p.get("value", 0.0)), dtype=float)

    if kind == "blocks":
        values = [float(v) for v in p.get("values", [0.0])]
        duration = float(p.get("block_duration", 3600.0))
        if duration <= 0:
            raise ValueError("block_duration must be positive")
        idx = ((t - t[0]) // duration).astype(int) % len(values)
        return np.asarray(values, dtype=float)[idx]

    if kind == "harmonic":
        mean = float(p.get("mean", 0.0))
        amplitude = float(p.get("amplitude", 1.0))
        period = float(p.get("period", 12.42 * 3600))
        phase = float(p.get("phase", 0.0))
        if period <= 0:
            raise ValueError("period must be positive")
        return mean + amplitude * np.sin(2 * np.pi * (t - t[0]) / period + np.deg2rad(phase))

    if kind == "linear":
        start = float(p.get("start", 0.0))
        end = float(p.get("end", 0.0))
        frac = (t - t[0]) / max(t[-1] - t[0], 1e-12)
        return start + (end - start) * frac

    if kind == "rotational":
        start = float(p.get("start", 0.0))
        rate = float(p.get("rate", 10.0))     # deg/hour
        return np.mod(start + rate * (t - t[0]) / 3600.0, 360.0)

    raise ValueError(f"unknown profile type '{kind}'")


def wind(tstart, tstop, dt, speed_spec, direction_spec):
    """Columns [t, speed, direction] for wind.txt."""
    t = time_axis(tstart, tstop, dt)
    u = np.clip(profile(t, speed_spec), 0.0, None)
    d = np.mod(profile(t, direction_spec), 360.0)
    return np.column_stack([t, u, d])


def tide(tstart, tstop, dt, level_spec):
    """Columns [t, water level] for tide.txt."""
    t = time_axis(tstart, tstop, dt)
    eta = profile(t, level_spec)
    return np.column_stack([t, eta])


def waves(tstart, tstop, dt, hs_spec, tp_spec):
    """Columns [t, Hs, Tp] for wave.txt."""
    t = time_axis(tstart, tstop, dt)
    hs = np.clip(profile(t, hs_spec), 0.0, None)
    tp = np.clip(profile(t, tp_spec), 0.0, None)
    return np.column_stack([t, hs, tp])


def write_series(path, data):
    np.savetxt(path, data, fmt="%0.6g")
