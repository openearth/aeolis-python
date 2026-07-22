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
piecewise:
    segments  {segments: [{type, duration, ...params}, ...]}
        Each segment covers *duration* seconds; a segment without a
        duration runs to the end of the series. Segment types are
        constant / linear / harmonic; a linear segment without an
        explicit start continues from the previous segment's end value.
        When every segment has a duration the series is generated only
        over the summed duration: AeoLiS repeats a boundary-condition
        file cyclically when it is shorter than the simulation
        (interp_circular / interp_circular_nearest).
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

    if kind == "segments":
        return _segments(t, p.get("segments") or [])

    raise ValueError(f"unknown profile type '{kind}'")


def _segments(t, segments):
    """Piecewise series: consecutive constant/linear/harmonic segments."""
    if not segments:
        raise ValueError("no segments defined")
    out = np.full(t.shape, np.nan, dtype=float)
    start = float(t[0])
    prev_val = None
    for i, seg in enumerate(segments):
        duration = seg.get("duration")
        # a segment without a duration runs to the end of the series
        if duration in (None, "", 0):
            end = float(t[-1])
            mask = t >= start
        else:
            end = start + float(duration)
            mask = (t >= start) & (t < end)
        if mask.any():
            ts = t[mask]
            kind = seg.get("type", "constant")
            if kind == "constant":
                vals = np.full(ts.shape, float(seg.get("value", 0.0)))
            elif kind == "linear":
                v0 = seg.get("start")
                v0 = float(v0) if v0 not in (None, "") else \
                    (prev_val if prev_val is not None else 0.0)
                v1 = float(seg.get("end", v0))
                frac = (ts - start) / max(end - start, 1e-12)
                vals = v0 + (v1 - v0) * np.clip(frac, 0.0, 1.0)
            elif kind == "harmonic":
                mean = float(seg.get("mean", 0.0))
                amplitude = float(seg.get("amplitude", 1.0))
                period = float(seg.get("period", 12.42 * 3600))
                phase = float(seg.get("phase", 0.0))
                if period <= 0:
                    raise ValueError("period must be positive")
                vals = mean + amplitude * np.sin(
                    2 * np.pi * (ts - start) / period + np.deg2rad(phase))
            else:
                vals = profile(ts, seg)
            out[mask] = vals
            prev_val = float(vals[-1])
        start = end
        if start >= t[-1]:
            break
    # samples beyond the covered duration: when every segment has a
    # duration the pattern is cyclic (AeoLiS repeats the file), so wrap
    # to the start of the cycle; otherwise hold the last value
    nan = np.isnan(out)
    if nan.any():
        good = ~nan
        total = segments_duration({"type": "segments", "segments": segments})
        if total and total > 0 and good.any():
            wrapped = t[0] + np.mod(t[nan] - t[0], total)
            idx = np.searchsorted(t[good], wrapped, side="right") - 1
            out[nan] = out[good][np.clip(idx, 0, int(good.sum()) - 1)]
        else:
            idx = np.where(good)[0]
            out[nan] = out[idx[-1]] if idx.size else 0.0
    return out


def segments_duration(spec):
    """Total covered time of a piecewise 'segments' spec in seconds, or
    ``None`` when any segment has no duration (i.e. runs to the end)."""
    if not isinstance(spec, dict) or spec.get("type") != "segments":
        return None
    segments = spec.get("segments") or []
    if not segments:
        return None
    total = 0.0
    for seg in segments:
        duration = seg.get("duration")
        if duration in (None, "", 0):
            return None
        total += float(duration)
    return total


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
