# v0.6.27

## New Features

- New potential EoS: Kolafa-Nezbeda LJ potential (`KolafaNezbeda`)
- New ideal EoS: `GlennJL`, that uses the NASA-Glenn coefficients provided by the package [`Glenn.jl`](https://github.com/ProfLeao/Glenn.jl)
- Cubic solver: speed and stability improvements, the function `cubic_poly` (cubic polynomial in $Z$) was replaced by a routine that uses first a cubic polynomial in $\eta = \frac{b}{V}$ first, and then refines the result (with a polynomial in $V$ for liquid an a polynomial in $Z - 1$ for the gas branch) if the error between the objective pressure and the current pressure is above a threshold.
- Association: general mixing rules can now be implemented by using `Clapeyron.ijab_mix(f,Δijab,Kij)`/ `Clapeyron.ijab_mix!(f,Δijab,Kij)`
- Association: the following association mixing functions are now public API: `bondvol_mix`,`epsilon_assoc_mix`,`dufal_mix`, `ijab_zero_mix`. Inplace versions are also implemented ( `!bondvol_mix`,`!epsilon_assoc_mix`,`!dufal_mix`, `!ijab_zero_mix`).
- `recombine!`: speed improvements in the case of GC models and association models.
- Flash methods: speed and handling improvements to `Tproperty`/`Pproperty`.
- Flash methods: new utility functions: `is_active_phase(result::FlashResult,i)`, `each_active_phase_index(result::FlashResult)`, to filter inactive (zero fraction or invalid phases)
- Flash methods: new method, `numphase(result::FlashResult,true)`, that returns the number of active phases.
- Flash methods: `RRTPFlash` and `MichelsenTPFlash` now return the unconverged K value in an indirect way (via the composition of the inactive phase).
- Flash methods: `tp_flash2` is now marked as public.
- SAFTgammaMie: misc database improvements
- ogUNIFAC: misc database improvements
- Speed improvements in database loading
