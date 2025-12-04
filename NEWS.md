# v0.6.20

## New Features

- revamp to Fugacity-based bubble/dew solvers. Improved speed of successive substitution iterations and support for Activity models.
- Activity models: support for second-order Michelsen TP flash. in VLE and LLE equilibria
- Activity models: support for QP/QT flash
- New method: thermodynamic factor (`thermodynamic_factor`)
- New method: `eos_repr`, to create parseable julia code capable of reproducing a model
- `BACKSAFT`: support for multicomponent mixtures.
- `crit_pure`: improved convergence.
- CoolProp: support for CoolProp 7.2
- Association: removed small static solvers for sizes 2-5 due to compilation slowdown.

## Method deprecations

- `ActivityDewTemperature` and `ActivityBubbleTemperature` were removed, `FugBubbleTemperature` and `FugDewTemperature` are now the default for activity methods, with proper suppor for non-condensables/non-volatiles

## Bug fixes

- Convergence failure in Michelsen TP flash when equilibria = :unkwown and LLE was detected.
- Fixes on `MultiphaseTPFlash`
