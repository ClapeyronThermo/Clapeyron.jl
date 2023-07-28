# v0.5.0

## New Features
- rework of `@newmodel`, `@newmodelgc` and `newmodelsingle` macros. Now they also define the outer constructor. for a simple EoS that does not require transformation of parameters, you can now do:
```julia
@newmodel MyModel EoSModel MyModelParam
#define locations for your model, relative to the current database location
Clapeyron.default_locations(::Type{MyModel}) = ["models/mymodel"]
#define references
Clapeyron.default_references(::Type{MyModel}) = ["TODO"]
```
- new macro `@newmodelsingleton`, that defines "singleton" EoSModels.
- Rework of MultiParameter EoS. there are two new `EoSModels` that represent Empiric, Multiparameter EoS:
    - `SingleFluid`: for single component fluids
    - `MultiFluid`: for multicomponent fluids, with specific mixing and departure rules
    `SingleFluid` and `MultiFluid` constructors are capable of parsing `CoolProp` JSON single component files. furthermore, you can use the available `CoolProp` single fluid library by just importing `CoolProp` into the current working enviroment (`using CoolProp`).
- Cubic models that have an alpha dependent of the acentric factor, can be built by passing `acentricfactor` directly to `userlocations`, instead of `alpha_userlocations` (#188)
- New Function: `RGas(model)` and `Rgas()`, that gives the value of the gas constant used by the model. defaults to `Clapeyron.R̄ = 8.31446261815324`
- New model: `HelmAct`, to use multiparameter EoS + activity coefficient models as the departure.
- New model: `XiangDeiters`
- New model: `TholLJ` (used in `LJRef`)
- New model: `EmpiricIdeal`, the ideal part of a multiparameter model, to be used in conjuction with other EoS
- New model: `AlyLeeIdeal`
- New model: `CPLNGEstIdeal` (http://dx.doi.org/10.1016/j.jngse.2014.04.011)
- CSVs allow defining a `sep` keyword in the inline CSV options:
```
Clapeyron Database File
my parameters [csvtype = like, sep = ;]
species,Mw
1,3,5,7-CYCLOOCTATETRAENE;156.22368
```

## Breaking changes
-  `IAPWS95`, `PropaneRef`, `Ammonia2023` are now of type `SingleFluid{EmpiricAncillary}`.
- `GERG2008`, `EOS-LNG` are now of type `MultiFluid{EmpiricAncillary,AsymmetricMixing,EmpiricDeparture}`
- `@newmodel` macros don't require defining external constructors anymore.
- `SpecialComp` (used by `pharmaPCSAFT`) is now a `ClapeyronParam` instead of an `EoSModel`
