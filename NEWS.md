# v0.5.0

## New Features

- Rework of MultiParameter EoS. there are two new `EoSModels` that represent Empiric, Multiparameter EoS:
    - `SingleFluid`: for single component fluids
    - `MultiFluid`: for multicomponent fluids, with specific mixing and departure rules
    `SingleFluid` and `MultiFluid` constructors are capable of parsing `CoolProp` JSON single component files. furthermore, you can use the available `CoolProp` single fluid library by just importing `CoolProp` into the current working enviroment (`using CoolProp`).
- cubic models that have an alpha dependent of the acentric factor, can be built by passing `acentricfactor` directly to `userlocations`, instead of `alpha_userlocations` (#188)
- New Function: `RGas(model)` and `Rgas()`, that gives the value of the gas constant used by the model. defaults to `Clapeyron.R̄ = 8.31446261815324`
- New model: `HelmAct`, to use multiparameter EoS + activity coefficient models as the departure.
- New model: `XiangDeiters`
- New model: `TholLJ` (used in `LJRef`)
- New model: `EmpiricIdeal`, the ideal part of a multiparameter model, to be used in conjuction with other EoS
- New model: `AlyLeeIdeal`
- New model: `CPLNGEstIdeal` (http://dx.doi.org/10.1016/j.jngse.2014.04.011)
- New model: translated-and-consistent Peng-Robinson (`tcPR`)
- New model: translated-and-consistent Redlich-Kwong (`tcRK`)
- New model: consistent PR - Twu (`cPR`)
- New alpha function: Twu-88 (`Twu88Alpha`)
 - New alpha function: soave-2019 (`Soave2019Alpha`)
## Breaking changes
-  `IAPWS95`, `PropaneRef`, `Ammonia2023` are now of type `SingleFluid{EmpiricAncillary}`.
- `GERG2008`, `EOS-LNG` are now of type `MultiFluid{EmpiricAncillary,AsymmetricMixing,EmpiricDeparture}`
