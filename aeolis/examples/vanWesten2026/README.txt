AeoLiS Ecomorphodynamic Framework: Simulation Setups
This folder contains the numerical model setups used to demonstrate the new ecomorphodynamic framework for AeoLiS base.

These setups accompany the publication: A New Ecomorphodynamic Framework for Process-Based Modeling that Incorporates Vegetation Dynamics in Coastal Dune Evolution
The repository is organized into five main directories corresponding to the simulation sets described in Section 3 of the manuscript.

01_growth_* : Vegetation Growth and Spreading
Conceptual 2D simulations (10x10 m) isolating vegetation development processes without wind or morphodynamic updating.
- [01_growth_a]: Local growth dominated by vertical development.
- [01_growth_b]: Clonal expansion showing faster horizontal spread with reduced vertical growth.
- [01_growth_c]: Seedling dispersal facilitating long-distance colonization (stepping-stones).
- [01_growth_d]: Inter-species competition illustrating spatial exclusion between two species.

02_shear_* : Local and Wake Shear Reduction
Semi-1D transect simulations (7.0 m) evaluating the sensitivity of vegetation-induced shear reduction and morphodynamic bed level change.
- [02_shear_a1 to a3]: Varying tiller height (h_veg).
- [02_shear_b1 to b3]: Varying tiller density (N_t).
- [02_shear_c1 to c3]: Varying drag efficiency (\beta_veg).
- [02_shear_d1 to d3]: Varying the downwind wake recovery parameter (c_1).
- [02_shear_e1 to e3]: Varying the dynamic lifting coefficient (\alpha_lift).
- [02_shear_f1 to f3]: Varying the static grain bouncing coefficient (b).

03_zeta : Sediment Skimming and Trapping
A single semi-1D simulation (60 m) demonstrating how the bed-interaction factor (\zeta) is computed across five varying surface types: bare sand, a non-erodible layer, short vegetation, tall vegetation, and a downwind wake region.

04_foredune_* : One-Dimensional Foredune Development
1D simulations (250 m) demonstrating the morphological impact of the combined vegetation components over a 5-year period under real-world forcing. Odd letters (a, c, e, g) simulate a tall/fast-growing species; even letters (b, d, f, h) simulate a short/slower-growing species.
- [04_foredune_a/b]: Baseline approach (no sediment skimming).
- [04_foredune_c/d]: Sediment skimming activated.
- [04_foredune_e/f]: Clonal expansion enabled.
- [04_foredune_g/h]: Seed dispersal enabled.

05_blowouts_* : Real-World Blowout Simulation
Decadal 2D simulations (800x800 m) modeling the evolution of five artificial foredune notches at Nationaal Park Zuid-Kennemerland (NPZK), Netherlands.
- [05_blowouts_a]: Main demonstration case reproducing observed ecomorphodynamic behavior.
- [05_blowouts_b]: Stronger vegetation development rates, leading to stabilization and premature notch closure.
- [05_blowouts_c]: Weaker vegetation development rates, leading to prolonged mobility and morphodynamic activity.
- [05_blowouts_d]: Sediment skimming disabled (\zeta=1.0), resulting in highly localized deposition directly at the vegetation edges.