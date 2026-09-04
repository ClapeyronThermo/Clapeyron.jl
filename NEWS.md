# v0.6.28

## New Features

- `Unitful.jl` extension: support for missing bulk methods
- `Unitful.jl` extension: support for `PH`, `PS`, `VT`, `QT`, `QP`, `TS` basis
- `Unitful.jl` extension: support for all flashes
- Experimental: `Clapeyron.MassFractions` wrapper to indicate that the vector of compositions is a mass composition vector.

## Bug fixes

- fix in cubic solver
- fix in AD with `TProperty`/`Pproperty`
- removed extra `show` in `MichelsenTPFlash` when material balance failed.
- removed extra gibbs evaluation when checking bubble/dew consistency
