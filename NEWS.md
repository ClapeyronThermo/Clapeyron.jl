# v0.6.30

## New Features

- New model: Liquid State Electrolyte model (`LS`)
- New models: Fe-ANN-Mie and Fe-ANN(s)-Mie machine learning Mie potentials (`FeANNMie`)
- SAFT-VRQ-Mie: improved speed, added support for first order Feynman-Hibbs Perturbation by passing the `fh_order = :fh1` keyword. The inner solvers now use the implicit function theorem to propagate derivatives.
- saturation: improved codes for spinodal saturation and isotherm exploring machinery.
- critical point: new methods: `mechanical_crit_bracket`, that finds a bracket that bounds the critical point. `x0_crit_bisection` is an alternative, opt-in initial point for `x0_crit_pure` and `x0_crit_mix`, for models that cannot be splitted or have critical points that cannot be calculated easily as just functions of their parameters.
- New iterator functions: `comps_iterator(model)`,`groups_iterator(model | group)`,`groups_iterator(model | group ,i)`, `sites_iterator(model | site, i)`, `iions_iterator(model | Z)`, `ineutral(model | Z)`, that are the corresponding function equivalents to `@comps`, `@groups`, `@sites`, `@iions`, `@ineutral` iterator macros.
- New function: `component_charges`, that returns the charges of each component. neutral models do not define this methods, so explicit solvent electrolyte methods can be distinguished just by adding a method to this function
- `has_sites` now checks if the model has a definition of `Clapeyron.getsites` instead of searching by parameters.
- new method: `critical_salt_fraction` for finding the critical salt fraction at which LLE is possible
- `MSAID`: better initial points, solver now uses implicit function theorem to propagate derivatives
- cubics: spinodal initial point now only requires one evaluation instead of 4

## Bug fixes

- `RackettLiquid`: typo in mixture values
- `AbbottVirial`: typo in constants
- `TsonopoulosVirial`: typo in constants
- `KU`: fixed critical point evaluation
- `MSAID`: fixed evaluation
- `debye_length`: fixed derivative propagation on zero partials
