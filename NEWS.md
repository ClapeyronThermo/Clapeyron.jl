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

## Bug fixes
- `SAFTVRMie` was allocating excesively because of unbound type parameter. 
- typos in `pharmaPCSAFT`
- `SanchezLacombe` didn't set `k` correctly when passed as `userlocations`
- `CPA`, SAFT equation of state and other EoS that implement association,don't need to specify `bondvol` and `epsilon_assoc`, when using non-associating species.
- correct implementation of `lb_volume` for `CPPCSAFT`
- better implementation of `lb_volume` for `pharmaPCSAFT`
