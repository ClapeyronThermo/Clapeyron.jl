# v0.6.11

## New Features

- The new minimum supported julia version is v1.10.
- Support for `ForwardDiff` v1.0.
- the following EoS now support parametric parameters:
  - `SAFTgammaMie` (including `structSAFTgammaMie`)
  - `PCPSAFT` (including `HeterogcPCPSAFT` and `HomogcPCPSAFT`)
  - `CPPCSAFT`
  - `sPCSAFT` (including `gcsPCSAFT`)
  - `pharmaPCSAFT`

Given a model with parametric parameters, one can now build another model with a different number type using the function `Clapeyron.promote_model(::Type{T},model) where T <: Number`.

- `SAFTgammaMie`: mixed segment paramters are now stored in the model parameters instead of inside the groups.
- Faster `split_model`
- Faster parameter instantiation, as now the `sources` and `sourcescsv` can be optionally `nothing`,and there is less copying of vectors.
- New Function: `USCT_temperature`. `USCT_mix` was renamed to `USCT_pressure` (The alias is still available, but it could be removed in future Clapeyron versions.)
- Estimation Framework: initial support for gradient optimization (using `ForwardDiff`) with parametric models.

## Bug Fixes

- fixed `CPA` initialization with custom parameters.
- fixed `CPA` default locations.
- fixed `ePCSAFT` initialization with custom ideal models.
- general flash: support for pure supercritical states.
- fix bugs in noncondensable/nonvolatiles Fugacity solver.
