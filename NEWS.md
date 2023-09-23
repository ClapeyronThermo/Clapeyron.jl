# v0.5.3
- Databases were standarized according to CAS. almost all components present in Clapeyron.jl databases are present in `@DB/properties/identifiers.csv`.
- COSMOSAC-2002 (`COSMOSAC02`),COSMOSAC-2010 (`COSMOSAC10`) and COSMOSAC-dispersion (`COSMOSACdsp`) can now read files from the NIST database found at https://github.com/usnistgov/COSMOSAC . to use those parameters, pass the keyword `use_nist_database = true`
- New model: doubly association perturbation theory (`DAPT`)
- New model: PCSAFT with association dependent hard sphere diameter (`ADPCSAFT`)
- New model: translated-and-consistent Peng-Robinson (`tcPR`)
- New model: translated-and-consistent Peng-Robinson, with Wilson and gE-res mixing rule (`tcPRW`)
- New model: translated-and-consistent Redlich-Kwong (`tcRK`)
- New Mixing Rule: residual (excess) gibbs energy mixing rule (`gErRule`)
- New model: consistent PR - Twu (`cPR`)
- New alpha function: Twu-88 (`Twu88Alpha`)
- New alpha function: soave-2019 (`Soave2019Alpha`)
- `DETPFlash` now supports the `equilibrium = :lle` option, to only find liquid phases.

# v0.5.2

## New Features
- Support for solid models (`Clapeyron.sle_solubility`, `Clapeyron.slle_solubility` and `Clapeyron.eutectic_point`).

## Bug fixes
- `eltype(::SAFTVRMieParam)` is defined.

# v0.5.1

## New Features
- Experimental support models with diferent types of parameters (that aren't `Float64`). This allows, among other things, creating models that have uncertainty in their parameters, and track that uncertainty across property calculations. `PCSAFT` and `SAFTVRMie` support this. (uncertainty support via `Measurements.jl` + `ForwardDiffOverMeasurements.jl` for autodiff rules)
- Models built automatically via the `@newmodel`, `@newmodelgc` and `newmodelsingle` macros now allow to pass single components as a string (`PR("water")`). There is also more compatibility with Group Contribution models `PCSAFT(["water" => ["H2O"=>1]],idealmodel = WalkerIdeal)`.
- minor inprovement on `x0_sat_pure` when the model cannot provide a virial coefficient.

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
