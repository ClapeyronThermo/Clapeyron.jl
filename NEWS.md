# v0.6.30

## New Features

- New model: Liquid State Electrolyte model (`LS`)
- Corresponding states models: revised formulation allows doing corresponding states between models with different gas constants. in particular, volume calculations, critical points and saturation pressures are calculated in the reference model first and then scaled. The API of corresponding states models was polished, so other models can subtype `CorrespondingStatesModel`. 
- `ECS`: the model type is just now a generic wrapper over a reference model and a scaling model, instead of a `SPUNG` specific model.
- `SPUNG`: moved over to a scaling type: `SPUNGScaling`, added exact solutions for single-component, soave-type cubics.
- New scaling models `NoECSScaling` (neutral scaling), `LJScaling` (vdW 1-fluid mixing rules) and `ChainLJScaling` (vdW 1-fluid mixing rules with segment). note that the chain scaling rule does not provide any additional contributions apart from considering the segment as a part of the scaling factors.
- New model types: Potential models, simple models for reference fluids. defined in molecular units (`T_scale == 1`, `Rgas = 1`).
- New potential model: `Pohl2023Mie6` for a ν-6 mie fluid (ν between 11 and 13)
- New potential model: `Johnson1993LJ`: LJ reference fluid used in `softSAFT`
- New potential model: `Thol2016LJ`: Before, the code used `TholLJ()`, that function now returns the `SingleFluid` stored inside `Thol2016LJ`
- New potential model: `KolafaNezbedaLJ`: potential model used by `KolafaNezbeda`. `KolafaNezbeda` is now a corresponding states model.
- New potential model: `PeTSLJ`: potential model used by `PeTS`. `PeTS` is now a potential model.
- `LJRef` now allows changing the inner LJ model.
- `softSAFT` now allows changing the inner LJ model. `softSAFT2016` is just an instance of `softSAFT` using `Thol2016LJ`.
- Lee-Kesler-Plöcker: new mixing formulation for `LKPSJT`, where the correlation for Zc interpolates the critical compressibility factors of the single and reference fluid instead of using the original correlation. The new formulation was selected as the default; the old correlation can still be used by passing the `version = 1` argument.
- New potential model: Fe-ANN-Mie and Fe-ANN(s)-Mie machine learning Mie potentials (`FeANNMie`)
- SAFT-VRQ-Mie: improved speed, added support for first order Feynman-Hibbs Perturbation by passing the `fh_order = :fh1` keyword. The inner solvers now use the implicit function theorem to propagate derivatives.
- saturation: improved codes for spinodal saturation and isotherm exploring machinery.
- critical point: new methods: `mechanical_crit_bracket`, that finds a bracket that bounds the critical point. `x0_crit_bisection` is an alternative, opt-in initial point for `x0_crit_pure` and `x0_crit_mix`, for models that cannot be splitted or have critical points that cannot be calculated easily as just functions of their parameters.
- New iterator functions: `comps_iterator(model)`,`groups_iterator(model | group)`,`groups_iterator(model | group ,i)`, `sites_iterator(model | site, i)`, `iions_iterator(model | Z)`, `ineutral(model | Z)`, that are the corresponding function equivalents to `@comps`, `@groups`, `@sites`, `@iions`, `@ineutral` iterator macros.
- New function: `component_charges`, that returns the charges of each component. neutral models do not define this methods, so explicit solvent electrolyte methods can be distinguished just by adding a method to this function
- `has_sites` now checks if the model has a definition of `Clapeyron.getsites` instead of searching by parameters.
- new method: `critical_salt_fraction` for finding the critical salt fraction at which LLE is possible
- `MSAID`: better initial points, solver now uses implicit function theorem to propagate derivatives
- cubics: spinodal initial point now only requires one evaluation instead of 4
- `enthalpy_vap`: Now supports any model with an available implementation of saturation pressure.
- `acentric_factor`: Now supports any model with an available implementation of saturation pressure.

## Bug fixes

- `RackettLiquid`: typo in mixture values
- `AbbottVirial`: typo in constants
- `TsonopoulosVirial`: typo in constants
- `KU`: fixed critical point evaluation
- `MSAID`: fixed evaluation
- `debye_length`: fixed derivative propagation on zero partials
- `SingleFluid`: fixed scaling of `lb_volume` and calls to saturation pressure at low temperatures.
- `x0_volume_gas`: specific model gas constant was not being considered during the virial volume solving.
