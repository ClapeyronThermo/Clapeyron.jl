# v0.6.4

## New Features

- New model: SAFT-VR-Mie with Gross-Vrabec quatrupolar contribution (`SAFTVRMieGV`)
- New model: Co-Oriented Fluid Functional Equation for Electrostatic interactions (`COFFEE`)
- Better support for evaluation of model properties at V == Inf (ideal gas limit)
- New method: `adiabatic_index`, that calculates the ratio between the isobaric and isochoric heat capacities.
- new API: `has_fast_crit_pure`, to indicate that models can calculate their pure critical point quickly. saturation initial guesses use the result of this function to decide if and when to call the `crit_pure` routine.
- speed ups in some pressure routines
-
## Bug fixes

- `MultiFluid` and `SingleFluid` models did not use the correct gas constant.
- Fix mixing rule in `SAFTVRMie`.
- `VT_identify_phase` now returns `:unknown` for an unstable state input.
- Typos in `TProperty` for pure models.

# v0.6.3

## New Features

- New Activity model: Hard-constraint Neural Network for Consistent Activity Coefficient Prediction(`HANNA`). requires the loading of the auxiliary package `ClapeyronHANNA`.
- New ideal model: PPDS correlations (`PPDSIdeal`)
- New Function: `reference_chemical_potential`, to calculate the reference chemical potential used in activity coefficient calculations
- New Functions: `activity` and `aqueous_activity`. both functions use the ability to change reference chemical potentials.
- New Functions: `spinodal_pressure` and `spinodal_temperature`.
- New function: `split_model_binaries`, that returns a list of all binary combinations of an n-component model.
- New Saturation Method: `CritExtrapolation`, that calculates saturation volumes and pressures via extrapolation from the critical point.
- `lb_volume` now has a three-arg version: `lb_volume(model,T,z)`.
- `p_scale` is now defined in terms of `T_scale` and `lb_volume`.
- Stability improvements for bubble/dew calculations.
- Speed improvements for single and multicomponent equilibria.
- Association solver is now faster for small association matrices.
- New association mixing rule: Mie-15 (`:mie15`,`:dufal`)
- Michelsen TP-Flash: in case of valid K values but single phase rachford-rice, the procedure will assume bubble or dew point as a first iteration.
- Joback: new submodule: `JobackGC` that provides all available properties using the joback correlations.
- SAFT-VR-Mie: speed improvements for calculation of association strengths.
- Cubics: Better initial point for single component saturation calculations.
- `split_model` now works for `ClapeyronParam`,`Symbol`,`Number`,`AbstractString`,`Tuple`,`Missing` and `Nothing`. before those could only be splitted if inside an `EoSModel`.
- `StructGroupParam` is deprecated, `GroupParam` has all the functionality of `StructGroupParam`.

## Bug fixes
- `SAFTgammaMie` fixes.
- `SingleFluid` has improved initial points for liquid volume evaluation.
- miscelaneous database improvements.
- `second_virial_coefficient` for cubics was ignoring the translation.
- improvements to the initial point of `SingleFluid`.

# v0.6.2

## Bug fixes
- `SAFTgammaMie` fixes.

# v0.6.1

## Bug fixes
- `tpd` fixes.
- `MultiPhaseTPFlash` fixes.
- Association: all temporary storage is now initialized

# v0.6.0

## New Features
- New models: Electrolyte models are now supported! We have introduced the `ESElectrolyte` framework which will let users combine any electrostatic model (`DH`, `MSA` and `Born`) and relative static permittivity model with any of our supported equations of state. Due to this flexibility, we now support four existing SAFT-type electrolyte equations (with planned support for more):
  - `ePCSAFT`
  - `SAFTVREMie`
  - `eSAFTVRMie`
  - `SAFTgammaEMie`
- New method: Two new methods specific to electrolytes have been added: `mean_ionic_activity_coefficient` and `osmotic_coefficient`, along with their saturated variants.
- New method: `MultiPhaseTPFlash`, that solves multiphase,multicomponent TP-flash with automatic phase detection. this method is now the default when calling `tp_flash` with more than two components and helmholtz-based models.
- New method: `Tproperty(model,p,prob,z,property)` to calculate temperatures, given pressure and another property.
- New model: To model solubility of salts, `SolidKs` has been added in order to obtain the solubility of salts using the infinite-dilution approach as opposed to the pure-fluid approach using `SolidHfus`.
- additional method: `x0_volume_liquid(model,p,T,z)` and `x0_volume_solid(model,p,T,z)` can be overloaded to calculate liquid an solid volumes, using the pressure as information. They are defined as `x0_volume_liquid(model,p,T,z) = x0_volume_liquid(model,T,z)` and `x0_volume_solid(model,p,T,z) = x0_volume_solid(model,T,z)`
- tangent plane distance (`tpd`) calculations are now faster.
- `VT_diffusive_stability` now uses `eigmin` instead of the full eigen calculation.
- `isstable` now works on (P,T,z) space, for the (V,T,z) space, use `VT_isstable`. there are now (P,T,z) versions of each stability function.
- calculation of volumes,saturation pressures and critical points of CPA models now defaults to the inner cubic model when there is no association present.
- The default association implementation now uses a combination of accelerated successive substitution and newton optimization. While increasing allocations, the method is faster.
- the default `volume` implementation now uses implicit AD to support derivatives. instead of propagating derivative information through the iterative procedure. This allows workloads of the type: `ForwardDiff.derivative(_p -> property(model,_p,T,z,phase = :l,vol0 = v0),p)` to be efficiently calculated.
- `Clapeyron.tpd` code has been optimized. `tpd` has new keywords: `break_first`, that tries to return a negative tpd as early as possible, `lle` for only calculating TPD in liquid phases, `strategy`, that changes the search strategy between a K-value search (`:wilson`), a pure component search (`:pure`) or both strategies (`default`).
- `Clapeyron.tpd` now supports activity models (if the keyword `lle` is set to `true`)
- `MichelsenTPFlash` will now try LLE if the initial equilibrium type is not set and VLE fails.
- New EoS: modified Lee-Kesler-Plöcker with consistent parameters (`LKPmod`)
- New EoS: Lee-Kesler-Plöker-equation of state, Sabozin-Jäger-Thol enhancement (`LKPSJT`, `enhancedLKP`)
## Bug fixes
- PCPSAFT: typo in unlike asssociation parameters

# v0.5.11

## New Features
- Support for reference states. A reference state is a point in V-T space where H = H₀ and S = S₀. Setting those have uses in Reaction equilibria and when comparing between different models.In particular, Reference states are stored in a `ReferenceState <: ClapeyronParam` in the ideal model parameters. The `BasicIdeal` model is, intentionally, the only ideal model in Clapeyron that does not have this struct and, as a consequence, it is not able to set reference states.
- Support for superancillaries via [`EoSSuperancillaries.jl`](https://github.com/ClapeyronThermo/EoSSuperancillaries.jl). When the package is loaded, initial saturation points for cubics and PCSAFT are overloaded to use superancillary evaluations instead of the general `x0_sat_pure` function. in the case of `PCSAFT` models, it also speeds up the evaluation of `crit_pure`.
- New EoS: EOS-CG (2021) (`EOS_CG`), a reference model for humid gases and CCS mixtures.
- New EoS: Lee-Kesler-Plöcker (`LKP`)
- New EoS: Shomate ideal model (`ShomateIdeal`)
- database: PCPSAFT,gcPCSAFT and gcPCPSAFT are updated to use the values of Rehner (2023).
- new functions: `helmholtz_free_energy_res`,`gibbs_free_energy_res`,`internal_energy_res`, `enthalpy_res`
- database: `ReidIdeal` now uses the poling coefficients by default.
- database: `JobackIdeal` has support for more common group fragments used in gcPCSAFT.
- `melting_temperature`, `sublimation_temperature` does not allocate anymore. Note that the function can still allocate if the EoS model itself allocates.

## bug fixes
- Incorrect value for CPA with water (#256)
- Bug in SAFT-VR-SW (#165)
- Bug in CP-PC-SAFT

# v0.5.10

## New Features
- Association models don't allocate anymore in the case of a single association site pair.
- `saturation_pressure(model,T)` (`ChemPotVSaturation,IsoFugacitySaturation`) does not allocate if the calculation does not require a critical point calculation. Note that the function can still allocate if the EoS model itself allocates. the same optimizations were applied to `saturation_temperature` (`AntoineSaturation`,`ClapeyronSaturation`), `sublimation_pressure` and `melting_pressure`.
- Bulk properties now accept a `vol0` initial point for the volume solver.
- SAFT-VR-Mie uses a divided approach for calculating `d`: if θ = ℂ*ϵᵢ/T > 1, then it uses a 10-point gauss-laguerre integrator. Otherwise, the Assen method of finding a cut point and integrating the rest is used. A description of the method is found here: https://teqp.readthedocs.io/en/latest/models/SAFT-VR-Mie.html. the cut allows for better accuracy at higher reduced temperatures.

## Bug fixes
- Peng-Robinson now uses more accurate `Ωa` and `Ωb` values
- CPA/sCPA now uses SI units as input.

# v0.5.9

## New Features
- New EoS: Solid SAFT-VR Mie (`SAFTVRSMie`)
- New EoS: Solid soft-SAFT (`solidsoftSAFT`)
- New property: sublimation pressure. `sublimation_pressure(model::CompositeModel,T)`
- New property: melting pressure. `melting_pressure(model::CompositeModel,T)`
- New property: sublimation temperature. `sublimation_temperature(model::CompositeModel,p)`
- New property: melting temperature. `melting_temperature(model::CompositeModel,p)`
- New property: triple point. `triple_point(model::CompositeModel)`
- `CompositeModel` was revamped to support more general equilibria. in particular it will be used to represent equilibria with Activity Models along with with Real Gases. As a result of these enhancements, `CompositeModel` now supports `bubble_pressure`,`bubble_temperature`,`dew_pressure`, and `dew_temperature`.
- `DETPFlash` supports LLE equilibria with activity models
- Bulk properties now accept a `vol0` initial point for the volume solver.

## Bug fixes
- `SAFTVRMie` was allocating excesively because of unbound type parameter.
- typos in `pharmaPCSAFT`
- `SanchezLacombe` didn't set `k` correctly when passed as `userlocations`
- `CPA`, SAFT equation of state and other EoS that implement association,don't need to specify `bondvol` and `epsilon_assoc`, when using non-associating species.
- correct implementation of `lb_volume` for `CPPCSAFT`
- better implementation of `lb_volume` for `pharmaPCSAFT`

# v0.5.8

## New Features
- `Base.getindex` and `Base.setindex` with `SingleParam`, `PairParam` and `AssocParam` now works with strings. the strings are compared with the components (or groups) stored in each param. in particular `AssocParam` allows set/get index methods if you pass a `Tuple{String,String}`:
```julia
julia> model = PPCSAFT(["water","ethanol"],assoc_options = AssocOptions(combining = :esd))
PPCSAFT{BasicIdeal} with 2 components:
 "water"
 "ethanol"
Contains parameters: Mw, segment, sigma, epsilon, dipole, dipole2, epsilon_assoc, bondvol

julia> model.params.bondvol[("water","a"),("water","b")]
0.35319

julia> model.params.bondvol[("water","a"),("water","b")] = 0.36
0.36

julia> model.params.bondvol[("water","a"),("water","b")]
0.36
```
- `PCPSAFT` is defined (alias for `PPCSAFT`)
- New EOS: Critical-point based PC-SAFT `CPPCSAFT` (https://doi.org/10.1021/ie502633e)
- Experimental: new `get_k`/`get_l`/`set_k!`/`set_l!` defined for cubics

## Bug Fixes
- bug in ether and aldehyde parameters in UNIFAC (https://github.com/ClapeyronThermo/Clapeyron.jl/issues/225)
- more strict checks for saturation temperature and better initial point.

# v0.5.7

## New Features
- New EoS: Polar PCSAFT with Quadrupolar interactions (`QPPCSAFT`)
- New EoS: Group-Contribution simplified PC-SAFT (`gcsPCSAFT`)
- New EoS: Group-Contribution homosegmented polar PC-SAFT (`gcPPCSAFT`)
- `Unitful.jl` support has been moved into an extension.
- `MultiComponentFlash.jl` support via extension: you can pass `model::EoSModel` to `MultiComponentFlash.flash_2ph`. there is also a new `MCFlashJL` method that calls `flash_2ph` using the `Clapeyron.tp_flash` interface.
- custom types can be passed to the `userlocations` keyword argument, defining `Clapeyron.can_nt(::datatype) = true` and `Clapeyron.to_nt(x::datatype)::Union{AbstractDict,NamedTuple}`

## Bug Fixes
- bug in `pharmaPCSAFT` mixing rules

# v0.5.6

## New Features
- `Clapeyron.diagvalues` now accepts `x::Number` (returning the same number)

## Bug Fixes
- more flexible sites parser (#214). Before, some site names where hardcoded.

# v0.5.5

## New Features
- Polar PC-SAFT now uses the Esper et al. (2023) parameters.

## Bug Fixes
- added promotion of types for solid solvers and `pharmaPCSAFT` (#212)

# v0.5.4

## New Features

- New Function: `export_model`, that allows exporting a current clapeyron model into a series of csv files.
- Metaheuristics.jl extension to allow passing a `Clapeyron.Estimator` to `Metaheuristics.optimize` (`optimize(f,estimator,method)`)
- `Clapeyron.Estimator` now has a `objective_form` field.

# v0.5.3

## New Features
- Databases were standardized according to CAS. almost all components present in Clapeyron.jl databases are present in `@DB/properties/identifiers.csv`.
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

# v0.4.13

## New Features

- (Experimental) Initial `Symbolics.jl` support for bulk properties on `EoSModel`. In particular, all Activity models, cubic models, SAFT models without association, and empiric models are supported. for example, this is now supported:

```julia
@variables y(..)[1:4], n(..)[1:4]
mixture = UNIFAC(["water","ethanol","methanol","1-propanol"])
moles = [n(t)[i] for i in 1:4] #Clapeyron accepts mole amounts, so it is not necessary to perform transformations to mole fractions
bc_l = y(t, 0.0) .~ activity_coefficient(mixture, 1.0, 298.15, moles)
```

## Bug Fixes

- automatic precompile is disabled in this version.

# v0.4.12

## New Features

- `@registermodel` macro is now a no-op. custom models that have a `components::Vector{String}` automatically support the printing interface and if they have `model.params.Mw`, they support molecular weight calculations.

- Cubic roots now use a real solver. this allows more stability and correctness in challenging EoS

# v0.4.11

## New Features

- The package now performs precompilation (via `PrecompileTools.jl`) of some commonly used functions and EoSModels. this will add more time at the installation of the package, in exchange with decreased loading times each time `Clapeyron` is loaded in a new julia session. you can turn on/off the feature with `Clapeyron.precompile_clapeyron!(setting::Bool)` (recomended in the case of developing the library). due to how precompilation is done, it is only done from julia 1.9 onwards.
- New EoS: NRTL with aspen paranetrization of τᵢⱼ : `aspenNRTL`
- `split_model` should be a little bit faster, and will perform correct splitting of group association models that use `get_group_idx`.

## Bug fixes

- Combining rule for epsilon in `SAFTgammaMie` updated from `√(ϵᵢ*ϵⱼ)*(σᵢ^3 * σⱼ^3)/σᵢⱼ^6` to `√(ϵᵢ*ϵⱼ*(σᵢ^3 * σⱼ^3))/σᵢⱼ^3`. the old behaviour can be obtained by passing the argument `epsilon_mixing::Symbol` to the model constructor. (#179)

# v0.4.10

## New Features

- `MichelsenTPFlash` now supports activity models, it also supports `CompositeModel` if they don't reach the multidimensional optimizer. with that, all combinations of 2-phase TP-Flash are supported in the following way:
    - Raoult: `CompositeModel`
    - Raoult with gas fugacity: `CompositeModel(components, gas = EoSModel)`
    - fugacity: any Helmholtz model
    - Activity + ideal gas: `Activity(components, puremodel = IdealModel)`
    - Activity + real gas: `Activity(components, puremodel = EosModel)` (`ActivityModel(components)` normally calls `ActivityModel(components,puremodel = PR)`)
- `RRTPFlash` now supports acceleration, non-condensables, non-volatiles, activity models and `CompositeModel`. (the same operations that `MichelsenTPFlash` supports.)
- `MichelsenTPFlash` and `RRTPFlash` provide initial guesses for LLE equilibria. `tp_flash(model,p,T,z,MichelsenTPFlash(equilibrium = :lle))` should suffice to calculate LLE flashes.
- `UNIFAC` models should be faster.

# v0.4.9

## New Features

- `ideal_consistency(model,V,T,z)` that checks if `da0/dV + sum(z)/V` is zero (or as close to zero as the Floating Point Format allows it.)
- broadcasting on `AssocParam` is defined (`bondvol .= 1e-10 .* bondvol .^3`)
- you can pass functions that build models instead of  EoSModel types. for example, now you can do:
    ```julia
    function myvdW(components;userlocations = String[],verbose = false)
        return vdW(components;userlocations = userlocations,verbose = verbose,alpha = SoaveAlpha)
    end

    model = Wilson(["water","ethanol"];puremodel=myvdW)
    ```

## Bug Fixes

- proper namespaces in `@registermodel` ([#161](https://github.com/ClapeyronThermo/Clapeyron.jl/issues/161))
- fixed bug in `MichelsenTPFlash` when using non-volatiles and `second_order = false`
- fixed bug when building `UNIFACFVPoly`

# v0.4.8

## New Features
- you can pass a named tuple to `userlocations` and its variants to build a model directly from those parameters (#156). for example, with PCSAFT:
```
julia> model = PCSAFT(["a1"],userlocations = (;
        Mw = [1.],
        epsilon = [2.],
        sigma = [3.],
        segment = [4.],
        k = [0.0;;], #matrix
        n_H = [1],
        n_e = [1],
        epsilon_assoc = Dict((("a1","e"),("a1","H")) => 1000.),
        bondvol = Dict((("a1","e"),("a1","H")) => 0.001)))
PCSAFT{BasicIdeal} with 1 component:
 "a1"
Contains parameters: Mw, segment, sigma, epsilon, epsilon_assoc, bondvol
```
- The `dense` option in `AssocOptions` is deprecated and it will be removed on 0.5.0. the sparse solver is not used anymore, the dense solver has performance advantages in all cases.
## Bug Fixes
- bug in `@registermodel` (#157)

# v0.4.7

## bug fixes
- `@newmodel`,`@newmodelsimple` and`@newmodelgc` macros had a bug where some structs and types didn't have the proper namespace (#154)

# v0.4.6

## Bug fixes
fixed bug on Michelsen TP flash introduced in the last release.

# v0.4.5

## Bug fixes
`DETPFlash` now returns results sorted by increasing molar volume (heavy liquid ->light liquid -> vapour) (fixes #151)

# v0.4.4

## New Features
- New EoS: Ammonia2023, an Empiric EoS with a new term (https://doi.org/10.1063/5.0128269)
- Speed improvements on some Empiric EoS
- Estimation Framework: the method is now specified as a csv option instead of just the title.
- Estimation Framework: support for setting species to be estimated.
## Bug fixes
- some references in the `SAFTVRMie`'s csv were incorrect.
- `recombine!` now works with `StructGroupParam`.

# v0.4.3

## Bug fixes
- incorrect indexing of sites when doing group-to-component transformation (e.g SAFTgammaMie) with more than one association site per molecule.
- better error handling when incorrect group types are assigned.

# v0.4.2

## Bug Fixes
- Typo in `assoc_translator` field when creating assoc sites with groups that have multiple association groups.
- Correct initialization of `β` in `rachfordrice`

# v0.4.1

## New Features

- New Clapeyron Parameter `StructGroupParam`, that stores groups, as well as it's bonds between them, used in heteronuclear GC approaches
- New EOS: structural SAFT-γ-Mie (`structSAFTgammaMie`)
- New EOS: heteronuclear GC-PC-SAFT (`gcPCSAFT`)
- in some limited cases, `split_model` can now differenciate between group and component parameters inside a model

# v0.4.0

## New Features

- New estimation framework, to fit new model parameters from known data. at the moment of this release, it should support all SAFTs, Cubics, Activity models and models that don't require any special pre-computed cache.

- `Base.transpose(model::EoSModel)` is defined. now this is valid code:
```julia
T = 300:350
comp = ["water"]
models = [m(comp) for m in (PR,PCSAFT,SAFTVRMie)]
saturation_pressure.(transpose(models),T) #produces a 51×3 Matrix{Tuple{Float64, Float64, Float64}}:
```

## Breaking Changes
- `x0_sat_pure` now returns `Vl` and `Vv` (in `[m^3]`) instead of `log10(Vl)` and `log10(Vv)`.
- in databases, `segment` is used instead of `m`
- in databases, `Vc` is used instead of `vc`
- in databases, `Pc` is used instead of `pc`
- in databases, `acentricfactor` is used instead of `w`


# v0.3.12

## new features
- new EoS: RKPR (`RKPR`)
- new functions: `cross_second_virial` and `equivol_cross_second_virial`, to calculate B12, at specified z conditions, or by setting equivolumetric mixing.
- cubic EoS now support `recombine!`
- `PenelouxTranslation` and `RackettTranslation` now cache their results, resulting in a speed up when using volume translated EoS.
- `Clapeyron.cite` now accepts the optional argument `out`, that can be `doi` (current default, DOI) or `:bib` (for BibTeX)
-
## Bug fixes
- various bug fixes in SAFT-VRQ-Mie. the Hard sphere term now performs the non-additive mixing rule for the multicomponent case. initializations with integer values are fixed.

# v0.3.11

## new features

- new EoS: PCSAFT with activity mixing rules (`GEPCSAFT(components; activity)`) (https://doi.org/10.1021/acs.iecr.2c03464)
- new EoS: UNIFAC-FV (`UNIFACFV`)
- new EoS UNIFAC-FV-poly (`UNIFACFVPoly`)
- `getparams` now supports inline csvs, custom locations and the ability to replace or swap out certain parameters, check `getparams` docs for more information
- new function: `recombine!` to recalculate combining and mixing rules after one parameter is modified. at the moment, SAFT and activity models have this defined.
- all quadratic mixing rules have an implace version (`sigma_LorentzBerthelot(sigma,zeta)` ->  `sigma_LorentzBerthelot!(sigma,zeta)`)
- `AssocParam` has `getindex`/`setindex!` methods defined.
- `GroupParam` has a new field, `grouptype`, to differenciate group parametrizations
-
## Deprecations

- `icomponents` was removed, use `Clapeyron.@comps` or `1:length(model)` instead
- `PairParam.diagvalues` was removed, use `Clapeyron.diagvalues(param)` instead

# v0.3.10

## New features
- `dew_pressure` and `dew_temperature` can now by calculated with `ActivityModel`s. options available can be passed via the new methods `ActivityBubblePressure`, `ActivityBubblwTemperature`, `ActivityDewPressure`, `ActivityDewTemperature`. helmholtz-based models also support solving `bubble_pressure` and `dew_pressure` using those methods (http://dx.doi.org/10.1021/acs.iecr.1c03800).
- There is support for γ-ϕ equilibria on `bubble_pressure` and `dew_pressure` with activity models, that is:, solving `yᵢϕᵢp = γᵢp₀ᵢϕ₀ᵢ𝒫ᵢ`
- New Correlation models: `LeeKeslerSat`, `DIPPR101Sat` for saturation pressure and temperature, `COSTALD`, `RackettLiquid`,`DIPPR105Liquid`, for saturated liquid volume.
- New models: Second Virial models (`AbbottVirial`,`TsonopoulosVirial` and `EoSVirial2`)
- New model: `CompositeModel`, used to hold saturation, liquid and vapor correlations. For example, instantiating an Activity Model, that uses Peng-Robinson for a gas model, DIPPR 101 for saturation and COSTALD for liquid volume, with a wilson activity coefficient, can be written as:
    ```
    julia> com = CompositeModel(["water","methanol"],liquid = COSTALD,saturation = DIPPR101Sat,gas = PR)
    Composite Model:
    Gas Model: PR{BasicIdeal, PRAlpha, NoTranslation, vdW1fRule}("water", "methanol")
    Liquid Model: COSTALD("water", "methanol")
    Saturation Model: DIPPR101Sat("water", "methanol")

    # Wilson activity model, using the composite model as the pure model
    julia> model = Wilson(["water","methanol"],puremodel = com)
    Wilson{CompositeModel{PR{BasicIdeal, PRAlpha, NoTranslation, vdW1fRule}, COSTALD, Nothing, DIPPR101Sat, Nothing}} with 2 components:
    "water"
    "methanol"
    Contains parameters: g, Tc, Pc, ZRA, Mw
    ```
- Improved initial points for `bubble_temperature` and `dew_temperature`.
- `MichelsenTPFlash` no uses an accelerated successive substitution method.

## Bug fixes

- cross-association weren't counted in some cases.
- incorrect initialization on `FugBubbleTemperature`
- `Clapeyron.lnϕ` uses the specialized algorithm for volume instead of the general one.

# v0.3.9

## New features
- New combining rules that are done at model creation time instead of runtime: `:cr1` and `:esd`/`elliott`. the old `elliott` combining rule was renamed to `elliott_runtime` (it is only used on `PharmaPCSAFT`, where is required)
## Bug fixes

- `crit_pure` on CPA was failing ([#112](https://github.com/ypaul21/Clapeyron.jl/issues/112), [#113](https://github.com/ypaul21/Clapeyron.jl/issues/113))
- error on bubble/dew pressure/temperature when using the `v0` starting point ([#113](https://github.com/ypaul21/Clapeyron.jl/issues/113))
- `crit_pure` was slower than expected on julia 1.8
- add error path on `SingleParam` when all data is missing

# v0.3.8

## New features

- `bubble_pressure`, `bubble_temperature`, `dew_pressure`, and `dew_temperature` can now support custom methods. (subtyping `BubblePointMethod` or `DewPointMethod`). the default methods are now named `ChemPotBubblePressure`, `ChemPotBubbleTemperature`,`ChemPotDewPressure`,`ChemPotDewTemperature`. you can provide some or all necessary initial conditions to those new methods:
    ```
    res = bubble_pressure(model,T,x,ChemPotBubblePressure(;y0;p0;vol0))
    ```
- New bubble/dew methods based on isofugacity conditions: `FugtBubblePressure`, `FugBubbleTemperature`,`FugDewPressure`,`FugtDewTemperature`.

- All bubble methods now support leaving some components out of the bubble phase,via the `nonvolatiles` keyword, supported by all available bubble solvers. similarly, all dew methods support the `noncondensables` keyword to leave some components out of the dew phase.

- `MichelsenTPFlash` can support both `nonvolatiles` and `noncondensables`.

- `SingleParam` and `PairParam` support indexing and broadcasting. `PairParam` has a modified indexing scheme, where `param[i]` will return the diagonal elements, `param[i,j]` will return off-diagonal elements. `setindex!` on `PairParam` is symmetric by default, that is, `param[i,j] = k` will also set `param[j,i]`. you can opt out of this by doing `param[i,j,false] = k`

- New function: `diagvalues`, for obtaining a vector view of the diagonal elements of a Pair or Single Parameter.

- New functions: `VT_gibbs_free_energy_res` and `VT_helmholtz_free_energy_res`

- New saturation method: `saturation_temperature(model,p,::SuperAncSaturation)`

- `index_reduction(model,x::AbstractVector{Bool})` will return a reduced model, according to the true values of the vector `x`. Passing another vector will perform the index reduction based on non-empty values, as usual.

- Initial (and optional) support for solid phase volume solving, via `x0_volume_solid`

- New EoS (experimental): `AnalyticalSLV`, that supports solid, liquid and vapour phases in one continous helmholtz functional.

- the database parser is faster and has better verbose reporting.

- `SAFTVRMie` and `SAFTGammaMie` should be a little faster on single component evaluations

## Bug Fixes

- `activity_coefficient` and `fugacity_coefficient` now works with `SVector` [#104](https://github.com/ypaul21/Clapeyron.jl/issues/104)
- there was a bug on calculating the length of GC models that went under `split_model`

## Deprecations

- `model.params.param.diagvalues` is not longer used in the codebase and it will be removed at the next breaking version. use `diagvalues(model.params.param)` or `model.params.param[i]` instead.

- `model.icomponents` is not longer used in the codebase and it will be removed at the next breaking version. `1:length(model)` or `Clapeyron.@comps` instead.

# v0.3.7

## New features
- You can now define custom saturation (temperature or pressure) solvers, via subtyping `SaturationMethod`. The default solver is now named `ChemPotVSaturation`.
- New methods for saturation pressure: `IsoFugacitySaturation`, `ChemPotDensitySaturation`, `SuperAncSaturation` (for cubics)
- All saturation methods support passing the `crit` keyword, to allow using precomputed critical points.
- All saturation methods that use a non linear solver support passing convergence options.
- The default solver for `saturation_temperature` was changed (now named `AntoineSaturation`). The new model is much faster and with a similar reliability. the old one is available via `ClapeyronSaturation`.
- New method for TP-Flash: `MichelsenTPFlash`
- New function: `tpd`, that tries to find all negative tangent plane distances
- New function: `lle_init`, that tries to find good LLE starting points, via `tpd`
- New function: `x0_saturation_temperature(model,p)`, that returns starting points for `saturation_temperature`.
- New function: `x0_psat(model,p)`, that returns starting points for `IsoFugacitySaturation`.
- New function: `antoine_coef` for reduced antoine coefficients. used to accelerate initial points, if available.
- New function: `cite(model)` returns all references used in the evaluation of a model.
- `bubble_temperature`, `dew_temperature`, `lle_temperature`, `azeotrope_temperature` and `vlle_temperature` methods now iterate directly over the extended-with-temperature non linear system that their pressure cointerparts.
- Cubics were revamped and generalized even more. those changes opened the way to introduce more non-standard cubics, all while being compatible with the initial guess procedures, advanced mixing rules and other Clapeyron optimizations for cubics.
- New cubic EoS: `Clausius`
- New cubic EoS: `Berthelot`
- New cubic EoS: Patel-Teja (`PatelTeja`)
- New cubic EoS: Patel-Teja-Valderrama (`PatelTeja`)
- New cubic EoS: Kumar-Upadhyay two parameter EoS (`KU`)
- New mixing rule: Modified Wong-Sandler (`ModWS`)
- `SoaveAlpha` now works with `vdW` and `PR` (returning the same value of `PRAlpha`)
- New EoS: Perturbed, truncated and shifted LJ EoS (`PeTS`)
- SPUNG is now a constructor for the new `ECS` model.
- `volume` now accepts a `vol0` initial point.
- `PCSAFT` should be a little bit faster for single component cases

# v0.3.6

## New Features

- `ZeroResidual` model (`a_res = 0`)

## Bug fixes

- `HVRule` uses the original rule instead of the modified one.
- `ismissingvalues` is preserved when a combining rule is applied. this allows us to track inputed parameters vs calculated ones.

In particular, this release is the one used in the `Clapeyron.jl` paper. and it's results should be reproducible by locking the version to v0.3.6.

# v0.3.5

## New Features
- new functions: `VT_mass_density`, `VT_mol_density`
- `vdW`, `RK` and `PR` models and variants use now `ABCubicParam`. `vdWParam`, `RKParam` and `PRParam` are now just aliases to  `ABCubicParam`.
- Roots.jl updated to 2.0
## Fixes
- fixes to UNIFAC models
- fixes to  `split_model` on groups

# v0.3.4

## New Features
- New EoS: Enhanced Predictive Peng Robinson 1978 (2022) (`EPPR78`)
- New EoS: Quantum Corrected Peng Robinson (2022) (`QCPR`)
- New EoS: Modified PCSAFT with T-dependent σ for water (passing `water08`) and k0 + k1T mixing rule (`pharmaPCSAFT`).
- New function: `index_reduction(model,z)` that removes components with compositions near zero. useful to reduce the dimensionality on multicomponent solvers.
- New function: `isstable(model,V,T,z)` to check for stability at the input conditions.
- New function: `gibbs_duhem(model,V,T,z)` to check if the Gibbs-Duhem relations holds for the input model at the specified conditions.
- New function: `kij_mix` and `pair_mix` for arbitrary pair combining rules
- volume solver can now perform aditional stability checks by passing `phase = :stable`, before only mechanical stability checks were done.
- association solver supports Elliot rule for combining rules, passed via `assoc_options(combining = :elliott`)
- association solver supports dense association matrix, passed via `assoc_options(combining = :dense_nocombining`)
- Faster `crit_pure`, `saturation_pressure`
- Faster `GERG2008` and `EOS_LNG`
- some cubic α-functions are faster for single component models
- improved initial saturation point for cubics
- `CPA` and `sCPA` now accept `AssocOptions`
- new mixing `SLKRule` for `SanchesLacombe`
- easier constructors for some empty models: `NoTranslation()` instead of `NoTranslation(NoTranslationParam())`

## Error fixes
- GC models can now be splitted in arbitrary component subsets
- `absolutetolerance` fields removed from almost all EoS (except COSMO models)
-  `second_virial_coefficient` now works with `SanchezLacombe`
- some starting points for `softSAFT` (original) and `BACKSAFT` were wrong. ``softSAFT2016` was not affected.
- `volume_compress` now accepts integer pressures, they are converted internally.
- `volume` for cubics now propagates the pressure type

# v0.3.3

## New Features
- `gibbs_solvation(model,T)` : For calculation of the solvation energy in a binary model
- speed up for activity coefficients and excess gibbs free energies for `NRTL`, `Wilson`, `UNIQUAC` and `UNIFAC`
- new EoS model: `softSAFT2016`, that uses the LJ helmholtz energy from thol et al. (2016)
- new utility function: `correct_composition_derivative`, that returns `∑xᵢμᵢ - G`. it should be zero or near zero for correctly written models.
## Error Fixes
- fix error in parameters for `softSAFT`
- fix errors in use of units. tests added.
- fix error on isobaric expansivity.

# v0.3.2

- fixes in volume solver
- fixes in initial guesses for `MultiFuidModel` and `CPA`
- all multicomponent solvers now accept an initial guess

# v0.3.1

Model fixes;
- `WSRule`
- `MonomerIdeal`
new functions:
`assoc_site_matrix(model::Union{SAFTModel,CPAModel},V,T,z)` to instantiate the association matrix
`ReidIdeal(model::JacobianIdeal)` to obtain a Reid ideal model from a Joback Ideal model
faster `crit_pure` on cubics

# v0.3.0

- Differential Evolution and Rachford-Rice Tp-flash algorithms have been implemented
- AssocOptions has been added to all SAFT equations
- Further speed-up to cubics
- sat_pure has been renamed saturation_pressure
- Numerous new multi-component VLE functions have been added (dew_pressure, bubble_temperature, VLLE_temperature etc.)
- New models: Full CPA, sCPA and ogUNIFAC

