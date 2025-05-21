# v0.6.12

## New Features

- SAFTgammaMie: easier constructor for inner SAFTVRMie model
- Clapeyron Parameters: support for typed empty constructor: (`SingleParam{BigFLoat}(name,components)`,`PairParam{BigFLoat}(name,components)`)
- XY Flash: added `QT`, `QP` and `VT` modules.
- XY Flash: support for second order properties with flash result if there is only one phase.
- NRTL: support for passing `tau` and `alpha` as input, instead of (`a`,`b`,`c`)
- CoolProp: initial support for superancillaries. At the moment, the superancillaries are used just as initial points. in future releases, We could return the result of the superancillary directly, to be in line with the CoolProp package.

## Bug Fixes

- Fixed incorrect value of `enthalpy_res`
- CoolProp: Support bigger buffer sizes.
- Implicit AD: misc bug fixes
- Fixed conversion of `MixedSegmentGCParam`
- Association: fix incorrect assumption of solved problem
- CPA: fixed initialization without `Pc`
- fix `promote_model` with `EoSVectorParam`
