# v0.6.18

## New Features

- `tp_flash`: improvements to `MichelsenTPFlash` in determination of the vapour phase fraction.
- `eutectic_point`: New keyword argument `x0` for passing an initial point.
- `eutectic_point`: support for `eutectic_point(model::SolidHfusModel,p)`
- `eutectic_point`: initial point calculation performs some succesive substitution iterations to improve the reliability of the solver (#466)
- support for `JutulDarcy` 0.3
