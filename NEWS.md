# v0.5.9

## New Features
- New EoS: Solid SAFT-VR Mie (`SAFTVRSMie`)
- New EoS: Solid soft-SAFT (`solidsoftSAFT`)
- New property: sublimation pressure. `sublimation_pressure(model::CompositeModel,T)`
- New property: melting pressure. `melting_pressure(model::CompositeModel,T)`
- New property: triple point. `triple_point(model::CompositeModel)`
- `CompositeModel` was revamped to support more general equilibria. in particular it will be used to represent equilibria with Activity Models along with with Real Gases. As a result of these enhancements, `CompositeModel` now supports `bubble_pressure`,`bubble_temperature`,`dew_pressure`, and `dew_temperature`.
- `DETPFlash` supports LLE equilibria with activity models
