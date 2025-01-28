# v0.6.8

## New Features

- New function: `partial_property`, for calculating partial properties at constant pressure and temperature.
- New functions: `widom_pressure`,`widom_temperature`,`ciic_pressure`,`ciic_temperature`, that calculate the maxima of isobaric heat capacity at constant pressure (widom) or constant temperature (CIIC).
- Ideal models: now all ideal models (with the exception of `BasicIdeal`) support optionally setting molecular weights.
- Database: the behaviour of `getparams(;return_sites)` was changed. Now, all single parameters used when building the `SiteParam` will be removed from the result (if sites are built.). The removed parameters will also not be checked for complete specification (non-existing sites are made equivalent to zero sites).
- Estimation: 
- ReidIdeal: Parameters for `a`,`b`,`c`,`d`,`e` are now included in `ReidIdealParam`.

## Bug Fixes

- Stability improvements in `xy_flash`.
- Fix K-value initialization when components are over JT temperature.
- Stability improvements for spinodal initialization, used for pure saturation pressure calculations.
