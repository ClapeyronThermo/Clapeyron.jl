# v0.6.14

## New Features

- Cubics: `ab_consts` is now defined automatically from `cubic_Δ`.
- Cubics: new abstract type for alpha models, `GeneralizedSoaveAlphaModel` (soave models with a different expression of `m` for each component.). `CPAAlphaModel`, `PTVAlpha` and `PatelTejaAlpha` are now subtypes of `GeneralizedSoaveAlphaModel`
- Electrolytes: revamp of code to make it more general.
- CoolProp: Experimental support for `CoolProp.PropsSI`
- new methods: `mass_entropy`,`mass_enthalpy`,`mass_internal_energy`,`mass_gibbs_free_energy`/`mass_gibbs_energy`, `mass_helmholtz_free_energy`/`mass_helmholtz_energy`, `mass_isochoric_heat_capacity`, `mass_isobaric_heat_capacity`. Note that those functions still take total volume and mol amounts as inputs.
- New alpha model: `LeiboviciAlpha`, a generalized soave model that works with any cubic (the value of `m` depends on the acentric factor and the values of `cubic_Δ`)
- `RKPR`: the `delta` parameter can now be specified from a database.
- `RKPRAlpha`: the `k1` and `k2` exponents can now be specified from a database
- `PatelTeja`: if `Vc` is not specified, it can be estimated from the acentric factor.
- `MTAlpha`: now works with `vdW`

## Bug Fixes

- typo in `QT` properties ([#380](https://github.com/ClapeyronThermo/Clapeyron.jl/issues/380))
- typo in `RKPRAlpha`
- improvements in spinodal calculation ([#382](https://github.com/ClapeyronThermo/Clapeyron.jl/discussions/382))
