# v0.4.10

## New Features

- `MichelsenTPFlash` now supports activity models, it also supports `CompositeModel` if they don't reach the multidimensional optimizer. with that, all combinations of 2-phase TP-Flash are supported in the following way:
    - Raoult: `CompositeModel`
    - Raoult with gas fugacity: `CompositeModel(components, gas = EoSModel)`
    - fugacity: any Helmholtz model
    - Activity + ideal gas: `Activity(components, puremodel = IdealModel)`
    - Activity + real gas: `Activity(components, puremodel = EosModel)` (`ActivityModel(components)` normally calls `ActivityModel(components,puremodel = PR)`)
- `RRTPFlash` now supports acceleration, non-condensables, non-volatiles, activity models and `CompositeModel`. (the same operations that `MichelsenTPFlash` supports.)
- `UNIFAC` models should be faster.
