# v0.6.15

## New Features

- Cubics: `cubic_ΔT` returns complex numbers, this allows evaluating cubics with non-factorizable atraction polynomials.
- New EoS: Yang-Frocsher-Richter Cubic (`YFR`)
- Symbolics.jl extension: support for residual properties
- Speed improvements in `Tproperty`, `Pproperty`,and chemical-potential based bubble/dew algorithms
- Speed improvements in XY-flashes for single-component models.

## Bug Fixes

- misc improvements to the database (removed duplicates and separate enantiomers)
- MultiComponentFlash.jl extension: fixed case of single phase result
