# v0.6.19

## New Features

- Activity models: support for second-order Michelsen TP flash. in VLE and LLE equilibria
- Activity models: support for tpd in VLE and LLE equilibria
- Activity models: new intrinsic: `lnγ_impl!(out,model,p,T,z)` that allows evaluation of activity coefficients without allocations
- Activity model performance improvements due to caching.
- New model: `EmpiricPseudoPure`: a Clapeyron implementation of CoolProp's pseudo pure models.
- New method: `RRQXFlash` for `qp_flash` and `qt_flash` first order (only used fugacity coeffients) secant roodfinding.
- `Clapeyron.tpd`: added ideal gas testing composition.

## Bug fixes

- Convergence failure in Michelsen TP flash when equilibria = :unkwown and LLE was detected.
- Fixes on `MultiphaseTPFlash`
