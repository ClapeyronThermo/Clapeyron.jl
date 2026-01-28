# v0.6.21

## New Features

- `MultiPhaseTPFlash`: support for Gamma-Phi models
- Cubics: translation now does not allocate in most cases.
- New model: translated industrial PC-SAFT (`iPCSAFT`)
- New model: VTPR with temperature-dependent translation (`TVTPR`)
- Electrolyte models: `mean_ionic_activity_coefficient` and `osmotic_coefficient` (and their saturated variants), can now be called directly without specifying salts, if the model forms binary salts. For example, (`osmotic_coefficient(model,p,T,m)`) is equivalent to(`osmotic_coefficient(model,salts,p,T,m)`), where `salts = Clapeyron.auto_binary_salts(model)`
- `DHModel`: improved numerical stability of `a_dh` when the electrolyte concentrations are small.
- association options can now be initialized directly from a symbol: (`PCSAFT(["water","ethanol"],assoc_options = :cr1)`)

## package deprecations

- JSON3.jl was removed in favour of JSON.jl

## Bug fixes

- `AntoineSaturation` now uses SI units
- fixes in `PenelouxTranslation` (inverted parameters for PR and RK models)
- Patel-Teja: the correlation for ζc is used first instead of the experimental critical volume
- Various documentation improvements
