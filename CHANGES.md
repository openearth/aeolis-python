# What changed on `aeolis-claude`

This branch starts from `bedinteraction_dev` at
[`4e5296d`](https://github.com/openearth/aeolis-python/commit/4e5296d) and changes only the model code
(`aeolis/`, `pyproject.toml`). All changes in one view:
[bedinteraction_dev...aeolis-claude](https://github.com/openearth/aeolis-python/compare/bedinteraction_dev...aeolis-claude).
The branch is for review, not for merging.

There are two kinds of commits:

- **Exact refactors.** They give results bitwise identical to `4e5296d`, checked on 13 test scenarios
  (blowout, foredune, barchan, parabolic, embryo dunes, Boschplaat, the regression cases and others).
- **`[physics]` commits.** They change results on purpose. Each one adds a configuration key, so the old
  behaviour stays available.

The default shear is still `method_shear = fft`. The new shear (`fft_grid`) has to be switched on.

## Where the gains come from

Wall time of whole runs on one machine (a 4-vcpu Xeon workstation), from `4e5296d` to this branch:

| Case | `4e5296d` | exact refactors only | + new shear, solver, time step | Speed-up |
|---|---:|---:|---:|---:|
| Blowout, 2 days (396 × 396) | 287 s | 160 s | 76 s | 3.8× |
| Foredune transect, 7 days | 42 s | 19 s | 2 s | 20× |
| Barchan, 7 days | 147 s | 137 s | 25 s | 5.9× |
| Parabolic dune, 7 days | 60 s | 42 s | 13 s | 4.6× |
| Embryo dunes, 7 days | 452 s | 360 s | 25 s | 18× |
| Boschplaat, 30 days (121 × 121) | 176 s | 95 s | 53 s | 3.3× |

The Boschplaat decade (ten years at hourly steps, one HPC core) takes about 5 h with `4e5296d`. The same
run takes 3.1–3.5 h with the exact refactors and the old shear, and 1.9–2.1 h with the new shear.

Ranked by effect:

1. **The shear on the model grid** (section 1) is the largest gain: 2–8× on its own. It is also the
   largest change in results.
2. **The bed-layer kernels** (section 2) are the largest exact gain. They give most of the factor 1.8
   between the first two columns.
3. **The steady-state solver** (section 3) adds 1.2–1.3× on large 2D grids. On the Boschplaat grid it
   saves almost nothing.
4. **The fused transport kernels** (section 4): transport and zeta become 3–6× faster per step.
5. **Start-up and output (netCDF)** (section 6) matter only for small and 1D runs, where they were most of
   the run time. For example, the 1D regression case went from 12 s to 2 s.
6. **The time-step limiter** (section 5) is a correctness fix, not a speed-up. Step counts stay about the
   same.

## 1. Shear perturbation on the model grid (`method_shear = fft_grid`) [physics]

The old `fft` method interpolates the bed onto a rotated, periodic grid, solves there and interpolates
back. That periodic box has a side effect at the edges. On a beach that rises 6–13 m across the domain,
the bed has to return to its starting height at the wrap, which acts as a cliff of that height. The cliff
gives a stress spike at the downwind edge and a wake upwind of it.

The new method solves the shear on the model grid itself. It rotates the wavenumbers instead of the grid,
so it needs no interpolation. Open boundaries are padded with the edge values. On the Boschplaat bed, the
rms error in τ/τ0 against a reference on a much larger domain drops from 0.19 (old) to 0.005 (new).

| Commit | What |
|---|---|
| [`17a99aa`](https://github.com/openearth/aeolis-python/commit/17a99aa) | The new method (first called `spectral`) |
| [`b8041cf`](https://github.com/openearth/aeolis-python/commit/b8041cf) | Cache the wavenumber factors (exact) |
| [`8257fcf`](https://github.com/openearth/aeolis-python/commit/8257fcf) | Separation bubble on the model grid |
| [`5ec9f78`](https://github.com/openearth/aeolis-python/commit/5ec9f78), [`c6d1554`](https://github.com/openearth/aeolis-python/commit/c6d1554) | No best-fit plane removal. Removing it deleted the real speed-up over a beach slope, and the option is gone |
| [`663d687`](https://github.com/openearth/aeolis-python/commit/663d687) | No short-wave filter by default. The filter removed the only damping of two-cell bed waves, which caused stripes on the Boschplaat foredune |
| [`29fc9e7`](https://github.com/openearth/aeolis-python/commit/29fc9e7) | Renamed `spectral` → `fft_grid`; `spectral` stays as an alias |
| [`e1379ab`](https://github.com/openearth/aeolis-python/commit/e1379ab) | `shear_transfer`: `weng` (one L), `weng_mode` (L per wavelength), `fca2010` (Fourrière, Claudin & Andreotti 2010) |
| [`3b32ab6`](https://github.com/openearth/aeolis-python/commit/3b32ab6) | `shear_short_waves`: optional treatment of waves a few cells long (default: none) |
| [`c8c6962`](https://github.com/openearth/aeolis-python/commit/c8c6962) | `shear_scale_A`, `shear_scale_B`: factors on the in-phase (A) and upwind-shift (B) coefficients |

Files: [`aeolis/processes/shear/spectral.py`](aeolis/processes/shear/spectral.py),
[`aeolis/processes/shear/separation.py`](aeolis/processes/shear/separation.py),
[`aeolis/wind.py`](aeolis/wind.py), new keys in [`aeolis/constants.py`](aeolis/constants.py).

`fca2010` with A and B scaled down (A × 0.4–0.5, B × 0.15–0.45) is currently being calibrated against the
Boschplaat decade (2014–2024 lidar). The defaults stay as they are until that is done.

## 2. Bed layers in compiled kernels (exact)

`bed.update` and `bed.mixtoplayer` looped over the layers in numpy for every cell. They now run as one
compiled kernel per column. Cells without pickup skip the bookkeeping.

| Commit | What |
|---|---|
| [`15a5c8e`](https://github.com/openearth/aeolis-python/commit/15a5c8e) | Per-column kernels (plus `setuptools<81` for `bmi-python`) |
| [`8fd5f3b`](https://github.com/openearth/aeolis-python/commit/8fd5f3b) | Same summation order as numpy, so results stay bitwise identical |
| [`4eba5e0`](https://github.com/openearth/aeolis-python/commit/4eba5e0) | Wired into `bed.update` and `bed.mixtoplayer` |
| [`cc3db0a`](https://github.com/openearth/aeolis-python/commit/cc3db0a) | Skip cells without pickup |

Files: [`aeolis/processes/bed/layers.py`](aeolis/processes/bed/layers.py), [`aeolis/bed.py`](aeolis/bed.py).

## 3. Steady-state solver

| Commit | What |
|---|---|
| [`e481d5c`](https://github.com/openearth/aeolis-python/commit/e481d5c) | The sweep iteration runs in compiled code (exact) |
| [`f6348f7`](https://github.com/openearth/aeolis-python/commit/f6348f7) | Sweep accepts any non-circular boundary name |
| [`865b499`](https://github.com/openearth/aeolis-python/commit/865b499) | **[physics]** `solver_order = upwind` (new default) |

The old solver sweeps the grid per wind quadrant and repeats until nothing changes, which took 20–60
sweeps per step on the blowout. The new order visits cells in the order of their upwind dependencies,
so most cells are solved exactly in one visit. Small recirculation loops are iterated locally.

Both orders converge to the same solution (within 3e-12 when converged tightly). The differences against
the old run are the old solver's own convergence error. On the blowout, solver time per step goes from
0.61 to 0.18 s. `solver_order = quadrant` gives the old sweeps.

Files: [`aeolis/processes/solvers/steadystate.py`](aeolis/processes/solvers/steadystate.py),
[`aeolis/advection.py`](aeolis/advection.py).

## 4. Fused transport and zeta kernels (round-off level)

Grain speed, the Newton solve and the equilibrium concentration took about 100 full-grid numpy passes per
step, and the vegetation part of zeta took several more. They now run as one compiled pass per cell. The
results agree with numpy to about 1e-12. `fused_kernels = F` selects the numpy code.

Per step on the blowout: transport 0.28 → 0.05 s, zeta 0.07 → 0.02 s.

| Commit | What |
|---|---|
| [`4fdac35`](https://github.com/openearth/aeolis-python/commit/4fdac35) | Fused kernels for grain speed, equilibrium transport and zeta |
| [`cf81581`](https://github.com/openearth/aeolis-python/commit/cf81581) | `duran_full` grain speed in parallel over rows (exact) |
| [`b814e8a`](https://github.com/openearth/aeolis-python/commit/b814e8a) | Okin leeside reduction in parallel (exact) |

Files: [`aeolis/processes/transport/kernels.py`](aeolis/processes/transport/kernels.py),
[`aeolis/processes/transport/zeta_kernels.py`](aeolis/processes/transport/zeta_kernels.py),
[`aeolis/transport.py`](aeolis/transport.py), [`aeolis/zeta.py`](aeolis/zeta.py).

## 5. Time-step limiter [physics]

The old limiter had five defects:

- `bed.initialize` overwrote `max_bedlevel_change`, so the configured value never applied.
- It tested the signed maximum of the bed change, so pure erosion was never limited.
- The bed change included avalanching, which does not scale with dt, so dt could only ratchet down.
- dt jumped straight back to its full value.
- Output times drifted, and the record at `tstop` was never written.

The new default `dt_limiter = aeolian` uses only the aeolian bed change and lands steps exactly on output
times and on `tstop`. `dt_limiter = legacy` gives the old behaviour. Defaults: `dt_growth = 0` (no growth
cap), `tstop_inclusive = F` (the old loop ran one step past `tstop`).

| Commit | What |
|---|---|
| [`45e6f8d`](https://github.com/openearth/aeolis-python/commit/45e6f8d) | The repaired limiter |
| [`3b1d0a5`](https://github.com/openearth/aeolis-python/commit/3b1d0a5) | `dt_growth` off by default |

Files: [`aeolis/model.py`](aeolis/model.py), [`aeolis/bed.py`](aeolis/bed.py), [`aeolis/constants.py`](aeolis/constants.py).

## 6. Start-up, output and threads (exact)

| Commit | What |
|---|---|
| [`a21c105`](https://github.com/openearth/aeolis-python/commit/a21c105) | Every numba kernel is cached; one type per scalar (an integer `dt` compiled a second version of every kernel); plotting imports only when needed (`import aeolis` 2.8 → 1.5 s) |
| [`3ae5283`](https://github.com/openearth/aeolis-python/commit/3ae5283) | netCDF output: one file open per write; no per-step copies for variables without statistics |
| [`d4bd5b2`](https://github.com/openearth/aeolis-python/commit/d4bd5b2) | Numba threads default to half the logical processors, with passive waiting (`AEOLIS_THREADS` overrides) |
| [`a76e51e`](https://github.com/openearth/aeolis-python/commit/a76e51e) | Skip the grid rotation when the grid angle is zero |

Files: [`aeolis/netcdf.py`](aeolis/netcdf.py), [`aeolis/model.py`](aeolis/model.py),
[`aeolis/__init__.py`](aeolis/__init__.py), and the kernel modules touched by `a21c105`.

## 7. Bug fixes and cleanup

| Commit | What |
|---|---|
| [`acab61a`](https://github.com/openearth/aeolis-python/commit/acab61a) | **[physics]** `euler_backward` silently ran Crank–Nicolson; it is now fully implicit |
| [`2260bc3`](https://github.com/openearth/aeolis-python/commit/2260bc3) | Legacy solver names (`trunk`, `pieter`, …) crashed; they now map to the current solvers with a warning |
| [`48e332c`](https://github.com/openearth/aeolis-python/commit/48e332c) | `1Dstacks` shear crashed; restart files were written in text mode; zero-step runs crashed; `tveg` missing in BMI |
| [`2938ab3`](https://github.com/openearth/aeolis-python/commit/2938ab3) | Dune erosion `Cs` had no default (KeyError) |
| [`d784341`](https://github.com/openearth/aeolis-python/commit/d784341) | Unused and non-working code removed (`advection_old.py` among others) |
| [`7733c91`](https://github.com/openearth/aeolis-python/commit/7733c91) | Regression references regenerated for the new defaults |

Found but **not fixed**, because each fix changes results:

- ds and dn are swapped in `gridparams` (matters only for rectangular cells).
- `np.max` is used where `np.maximum` is meant in `erosion.py`.
- z0 is indexed with a Fourier mask in `shear.py` (matters only for z0 that varies in space).
- `Cu` and `Ct` are scaled by `accfac`.
- The separation bubble depends on the grid spacing.

## Not on this branch

The benchmark harness, HPC scripts, validation cases and new unit tests are kept separately. Some commit
messages refer to them, for example the "gate" that checks bitwise equality with `4e5296d`.
