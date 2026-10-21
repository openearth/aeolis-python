# Overview of aeolis-claude

This branch makes AeoLiS faster and fixes a number of errors in the model code. It starts from the
`bedinteraction_dev` branch (commit `4e5296d`) and changes only the model itself: the `aeolis/` folder
and `pyproject.toml`. The branch is meant for review and is not meant to be merged.

[This comparison](https://github.com/openearth/aeolis-python/compare/bedinteraction_dev...aeolis-claude)
shows every changed line. This document explains which changes matter and why.

## Two kinds of changes

Most commits make the code faster and leave the results unchanged. We ran 13 test cases with both
versions and compared every output variable. The output is identical to the last bit. The test cases
include a blowout, a foredune, a barchan, a parabolic dune, embryo dunes, the Boschplaat on Terschelling
and the regression cases of the test suite.

The other commits change the results on purpose. They fix an error or add a better method. Their commit
messages start with `[physics]`. Each of them comes with a configuration key, so you can still switch
back to the old behaviour.

## Improvements

The table shows the run time of six cases on one workstation. The first column is
`bedinteraction_dev`. The second column adds only the changes that keep the results the same. The third
column also includes the changes that alter results: the new shear method (`method_shear = fft_grid`),
the new solver order and the repaired time step.

| Case | `bedinteraction_dev` | Same results | All changes | Speed-up |
|---|---:|---:|---:|---:|
| Blowout, 2 days (396 × 396 cells) | 287 s | 160 s | 76 s | 3.8× |
| Foredune transect, 7 days | 42 s | 19 s | 2 s | 20× |
| Barchan, 7 days | 147 s | 137 s | 25 s | 5.9× |
| Parabolic dune, 7 days | 60 s | 42 s | 13 s | 4.6× |
| Embryo dunes, 7 days | 452 s | 360 s | 25 s | 18× |
| Boschplaat, 30 days (121 × 121 cells) | 176 s | 95 s | 53 s | 3.3× |

For a longer run: ten years of the Boschplaat at hourly steps take about 5 hours on one core of the
Deltares cluster with `bedinteraction_dev`. The same run takes 3.1 to 3.5 hours with the changes that
keep the results the same, and about 2 hours with the new shear method.

### Which changes matter most

1. **Wind shear on the model grid** (section 1) gives the largest gain. It makes the model 2 to 8 times
   faster on its own. It also changes the results more than any other commit.
2. **Bed layers in compiled code** (section 2) give the largest gain without any change in results.
3. **The new solver order** (section 3) saves another 20 to 30 % on large 2D grids. On the 121 × 121
   Boschplaat grid it saves almost nothing.
4. **Combined transport calculations** (section 4) make the transport and zeta steps 3 to 6 times faster.
5. **Start-up and output** (section 6) matter only for small and 1D runs. There, they took most of the
   run time. The 1D regression case now runs in 2 s instead of 12 s.
6. **The repaired time step** (section 5) fixes errors. It hardly changes the run time.

## 1. Wind shear on the model grid

Dunes change the wind. The wind speeds up over a crest and slows down in front of a dune. AeoLiS
computes this effect with the analytical solution of Weng et al. (1991). It splits the bed into waves
with a Fourier transform (FFT), computes the shear stress response to each wave and adds the responses
up.

The existing method (`method_shear = fft`) does this on a separate grid that is rotated with the wind.
The model copies the bed onto that grid, computes the shear there and copies the result back. This has
two drawbacks:

- The copying between the grids smooths and shifts the shear field.
- The rotated grid wraps around: its downwind edge connects to its upwind edge. On a beach that rises
  6 to 13 m across the domain, this connection acts as a cliff. The cliff produces a false peak in shear
  stress at the downwind edge and a false sheltered zone upwind of it.

The new method (`method_shear = fft_grid`) computes the shear directly on the model grid. It rotates the
wave directions in the Fourier transform instead of rotating the grid. Beyond the open edges, it extends
the bed with the edge values, so no cliff appears. We compared both methods with the same calculation on
a much larger domain, which serves as the reference. On the Boschplaat bed, the error in shear stress
drops from 19 % with the old method to 0.5 % with the new one.

The new method is optional. The default stays `method_shear = fft`. We are calibrating the new method
against the Boschplaat lidar data of 2014 to 2024 before we propose it as the default.

| Commit | What it does |
|---|---|
| [`17a99aa`](https://github.com/openearth/aeolis-python/commit/17a99aa) | Adds the new method. Its first name was `spectral`. |
| [`8257fcf`](https://github.com/openearth/aeolis-python/commit/8257fcf) | Computes the flow separation behind steep lee slopes on the model grid too, and supports periodic boundaries. |
| [`6bdc8fe`](https://github.com/openearth/aeolis-python/commit/6bdc8fe), [`5ec9f78`](https://github.com/openearth/aeolis-python/commit/5ec9f78), [`c6d1554`](https://github.com/openearth/aeolis-python/commit/c6d1554) | The first version subtracted the average slope of the bed before the Fourier transform. That removed the real speed-up of the wind over a sloping beach. The later two commits switch it off and remove the option. |
| [`663d687`](https://github.com/openearth/aeolis-python/commit/663d687) | Switches off the smoothing filter by default. The filter removed the response to bed waves of two to three cells. Those waves then no longer flattened out and grew into stripes on the Boschplaat foredune. |
| [`29fc9e7`](https://github.com/openearth/aeolis-python/commit/29fc9e7) | Renames `spectral` to `fft_grid`. The old name still works. |
| [`e1379ab`](https://github.com/openearth/aeolis-python/commit/e1379ab) | Adds `shear_transfer`, the formula for the response to each wave. `weng` is the existing formula with one length scale L. `weng_mode` uses the wavelength itself as L, so L is no longer a tuning parameter. `fca2010` uses the fit of Fourrière, Claudin and Andreotti (2010). |
| [`3b32ab6`](https://github.com/openearth/aeolis-python/commit/3b32ab6) | Adds `shear_short_waves`, an optional damping of waves only a few cells long. The default is `none`. |
| [`c8c6962`](https://github.com/openearth/aeolis-python/commit/c8c6962) | Adds `shear_scale_A` and `shear_scale_B`. A sets how much the shear stress rises on a crest. B sets how far upwind of the crest the stress peaks; this shift is what makes dunes grow. The keys multiply A and B of the chosen formula. |
| [`b8041cf`](https://github.com/openearth/aeolis-python/commit/b8041cf) | Speeds up the existing `fft` method. It computes the parts that depend only on the grid once instead of every step. Same results. |

Main files: [`aeolis/processes/shear/spectral.py`](aeolis/processes/shear/spectral.py),
[`aeolis/processes/shear/separation.py`](aeolis/processes/shear/separation.py),
[`aeolis/wind.py`](aeolis/wind.py) and [`aeolis/constants.py`](aeolis/constants.py).

## 2. Bed layers in compiled code

AeoLiS stores the bed as a stack of layers, each with its own grain size distribution. Every time step,
`bed.update` and `bed.mixtoplayer` move sediment between these layers. The old code did this with numpy
operations over the whole grid, layer by layer. The new code loops over each column of layers in
compiled code (numba). The results stay identical to the last bit.

| Commit | What it does |
|---|---|
| [`15a5c8e`](https://github.com/openearth/aeolis-python/commit/15a5c8e) | Adds the compiled functions. It also pins `setuptools<81`, which `bmi-python` needs. |
| [`8fd5f3b`](https://github.com/openearth/aeolis-python/commit/8fd5f3b) | Adds up numbers in the same order as numpy, so the results match to the last bit. |
| [`4eba5e0`](https://github.com/openearth/aeolis-python/commit/4eba5e0) | Uses the compiled functions in `bed.update` and `bed.mixtoplayer`. |
| [`cc3db0a`](https://github.com/openearth/aeolis-python/commit/cc3db0a) | Skips cells where no sand was picked up. |

Main files: [`aeolis/processes/bed/layers.py`](aeolis/processes/bed/layers.py) and [`aeolis/bed.py`](aeolis/bed.py).

## 3. Solver order

The steady-state solver computes the sediment concentration in the air. Each cell depends on its upwind
neighbours. The old solver sweeps the grid once for each wind quadrant and repeats the sweeps until the
concentrations stop changing. On the blowout, this takes 20 to 60 sweeps per time step.

The new solver visits the cells in upwind order: each cell comes after the cells it depends on. Most
cells are then solved in a single visit. Where the wind circles back, the solver iterates only over that
small group of cells.

Both solvers converge to the same solution. Small differences remain, because the old solver stops at
its convergence limit. On the blowout, the solver time per step drops from 0.61 s to 0.18 s. With
`solver_order = quadrant` you get the old sweeps back.

| Commit | What it does |
|---|---|
| [`e481d5c`](https://github.com/openearth/aeolis-python/commit/e481d5c) | Runs the old sweeps in compiled code. Same results. |
| [`f6348f7`](https://github.com/openearth/aeolis-python/commit/f6348f7) | Lets the sweeps accept every boundary type, as the rest of the model does. |
| [`865b499`](https://github.com/openearth/aeolis-python/commit/865b499) | Adds the new solver order and makes it the default (`solver_order = upwind`). |

Main files: [`aeolis/processes/solvers/steadystate.py`](aeolis/processes/solvers/steadystate.py) and
[`aeolis/advection.py`](aeolis/advection.py).

## 4. Combined transport calculations

Grain speed, equilibrium transport and the vegetation part of zeta used about 100 separate numpy
operations per time step. Each operation passed over the whole grid. Now one compiled loop does all of
it per cell. On the blowout, transport takes 0.05 s per step instead of 0.28 s, and zeta takes 0.02 s
instead of 0.07 s.

The results differ from the numpy version by about 1e-12. These differences come from rounding.
`fused_kernels = F` switches back to the numpy version.

| Commit | What it does |
|---|---|
| [`4fdac35`](https://github.com/openearth/aeolis-python/commit/4fdac35) | Combines grain speed, equilibrium transport and zeta into compiled loops. |
| [`cf81581`](https://github.com/openearth/aeolis-python/commit/cf81581) | Computes the `duran_full` grain speed for several grid rows at once. Same results. |
| [`b814e8a`](https://github.com/openearth/aeolis-python/commit/b814e8a) | Computes the Okin lee-side reduction for several grid rows at once. Same results. |

Main files: [`aeolis/processes/transport/kernels.py`](aeolis/processes/transport/kernels.py),
[`aeolis/processes/transport/zeta_kernels.py`](aeolis/processes/transport/zeta_kernels.py),
[`aeolis/transport.py`](aeolis/transport.py) and [`aeolis/zeta.py`](aeolis/zeta.py).

## 5. Time step

AeoLiS shortens the time step when the bed changes too fast (`max_bedlevel_change`). This limiter had
five errors:

- `bed.initialize` overwrote `max_bedlevel_change`, so the value from the configuration file never
  applied.
- The limiter looked at the largest bed rise only, so it never limited erosion.
- It included avalanching in the bed change. Avalanching does not scale with the time step, so the time
  step could only shrink.
- After a short step, the time step jumped straight back to its full length.
- The output times drifted, and the model never wrote the output at `tstop`.

The new limiter (`dt_limiter = aeolian`) looks only at the bed change from wind transport. It also
shortens steps so that they end exactly on every output time and on `tstop`. `dt_limiter = legacy`
gives the old behaviour.

| Commit | What it does |
|---|---|
| [`45e6f8d`](https://github.com/openearth/aeolis-python/commit/45e6f8d) | Adds the new limiter and makes it the default. The model no longer runs one step past `tstop` (`tstop_inclusive = F`). |
| [`3b1d0a5`](https://github.com/openearth/aeolis-python/commit/3b1d0a5) | Switches off the limit on how fast the time step may grow again (`dt_growth = 0`). |

Main files: [`aeolis/model.py`](aeolis/model.py), [`aeolis/bed.py`](aeolis/bed.py) and
[`aeolis/constants.py`](aeolis/constants.py).

## 6. Start-up and output

| Commit | What it does |
|---|---|
| [`a21c105`](https://github.com/openearth/aeolis-python/commit/a21c105) | Stores every compiled function on disk, so the next run does not compile it again. Passes `dt` always as a float; an integer `dt` used to compile a second version of every function. Loads the plotting libraries only when needed: `import aeolis` takes 1.5 s instead of 2.8 s. |
| [`3ae5283`](https://github.com/openearth/aeolis-python/commit/3ae5283) | Opens the netCDF file once per output time instead of twice, and stops copying output variables that have no statistics. |
| [`d4bd5b2`](https://github.com/openearth/aeolis-python/commit/d4bd5b2) | Uses half the logical processors by default and lets idle threads sleep instead of wait actively. The environment variable `AEOLIS_THREADS` sets another number. |
| [`a76e51e`](https://github.com/openearth/aeolis-python/commit/a76e51e) | Skips the grid rotation when the grid is not rotated. |

All of these keep the results the same. Main files: [`aeolis/netcdf.py`](aeolis/netcdf.py),
[`aeolis/model.py`](aeolis/model.py) and [`aeolis/__init__.py`](aeolis/__init__.py).

## 7. Bug fixes

| Commit | What it does |
|---|---|
| [`acab61a`](https://github.com/openearth/aeolis-python/commit/acab61a) | Makes the `euler_backward` solver fully implicit. It used to run Crank–Nicolson without saying so. This changes results. |
| [`2260bc3`](https://github.com/openearth/aeolis-python/commit/2260bc3) | Old solver names such as `trunk` and `pieter` crashed the model. They now map to the current solvers with a warning. |
| [`48e332c`](https://github.com/openearth/aeolis-python/commit/48e332c) | Fixes four crashes and errors: `1Dstacks` shear, restart files written as text, runs of zero steps and the missing `tveg` in the BMI interface. |
| [`2938ab3`](https://github.com/openearth/aeolis-python/commit/2938ab3) | Gives the dune erosion coefficient `Cs` a default value again. Without it, the model stopped with a KeyError. |
| [`d784341`](https://github.com/openearth/aeolis-python/commit/d784341) | Removes code that nothing used or that did not work, such as `advection_old.py`. |
| [`7733c91`](https://github.com/openearth/aeolis-python/commit/7733c91) | Updates the reference output of the regression tests to the new defaults. |

## Errors found but not fixed yet

Each of these fixes would change results, so they need a decision first:

- `gridparams` swaps ds and dn. This matters only for cells that are not square.
- `erosion.py` uses `np.max` where it should use `np.maximum`.
- `shear.py` indexes z0 with a Fourier mask. This matters only when z0 varies in space.
- `Cu` and `Ct` are scaled by `accfac`.
- The size of the separation bubble depends on the grid cell size.

## What is not on this branch

We developed this work on a separate branch that also holds the test cases, the comparison scripts, the
cluster scripts and new unit tests. This branch leaves them out, so that it shows only the model code.
Some commit messages still refer to them, for example to the "gate", the script that checks that results
stay identical to `bedinteraction_dev`.
