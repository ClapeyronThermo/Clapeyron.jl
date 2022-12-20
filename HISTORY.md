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