# v0.6.9

## New Features

- `pharmaPCSAFT` now considers `water` and `water08` for it's T-dependend hard-sphere diameter. Before, only `water08` used the specialized behaviour, and `water` returned the stock PCSAFT parameters.
- speed improvements in `MichelsenTPFlash` and fugacity-based bubble/dew calculations

## Bug Fixes

- Fix errors with the use of `EoSSuperancillaries` along with `pharmaPCSAFT`
- rachford-rice: Fix errors when there are NaN K-values corresponding to non-condensables or non-volatiles.
- `tp_flash`: Fix errors when there are NaN K-values corresponding to non-condensables or non-volatiles.
- typo in `x0_sat_pure_spinodal`
