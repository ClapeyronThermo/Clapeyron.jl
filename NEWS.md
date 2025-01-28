# v0.6.8

## New Features

- New function: `partial_property`, for calculating partial properties at constant pressure and temperature.
- New functions: `widom_pressure`,`widom_temperature`,`ciic_pressure`,`ciic_temperature`, that calculate the maxima of isobaric heat capacity at constant pressure (widom) or constant temperature (CIIC).
- New models: original UNIFAC 2.0 (`ogUNIFAC2`, from [doi.org/10.1016/j.cej.2024.158667](https://doi.org/10.1016/j.cej.2024.158667)) and modified (Dortmund) UNIFAC 2.0 (`UNIFAC2`, from [doi.org/10.48550/arXiv.2412.12962](https://doi.org/10.48550/arXiv.2412.12962)).
- Ideal models: now all ideal models (with the exception of `BasicIdeal`) support optionally setting molecular weights.
- Database: the behaviour of `getparams(;return_sites)` was changed. Now, all single parameters used when building the `SiteParam` will be removed from the result (if sites are built.). The removed parameters will also not be checked for complete specification (non-existing sites are made equivalent to zero sites).
- Estimation: more flexibility in setting indices.
- ReidIdeal: Parameters for `a`,`b`,`c`,`d`,`e` are now included in `ReidIdealParam`.

## Bug Fixes

- Stability improvements in `xy_flash`.
- Fix K-value initialization when components are over JT temperature.
- Stability improvements for spinodal initialization, used for pure saturation pressure calculations.
