# v0.6.20

## New Features

- revamp to Fugacity-based bubble/dew solvers. Improved speed of successive substitution iterations and support for Activity models.
- Activity models: support for second-order Michelsen TP flash. in VLE and LLE equilibria
- Activity models: support for QP/QT flash
- Improved implicit differentiation for all solvers.
- New method: thermodynamic factor (`thermodynamic_factor`)
- New method: `eos_repr`, to create parseable julia code capable of reproducing a model
- `BACKSAFT`: support for multicomponent mixtures.
- `crit_pure`: improved convergence.
- CoolProp: support for CoolProp 7.2
- Association: removed small static solvers for sizes 2-5 due to compilation slowdown.
- Cubics: improved robustness in the volume solver
- PCSAFT: improved robustness in the volume solver at low temperatures

## Method deprecations

- `ActivityDewTemperature` and `ActivityBubbleTemperature` were removed, `FugBubbleTemperature` and `FugDewTemperature` are now the default for activity methods, with proper support for non-condensables/non-volatiles
- `DETPFlash` now uses another global optimizer, self-adaptive spherical search algorithm, with the main intention of dropping `BlackBoxOptim` as a dependency.

## Bug fixes

- Convergence failure in Michelsen TP flash when `equilibrium = :unknown` and LLE was detected.
- Fixes on `MultiphaseTPFlash`
- Various fixes on `saturation_pressure` initial points
- Fixed an extra `RT` division in `Obj_de_tp_flash`; it now returns the molar reduced Gibbs energy `G/(nRT)`.
- `__eval_G_DETPFlash(model::EoSModel, p, T, ni, equilibrium)` now returns reduced Gibbs energy, consistent with the `PTFlashWrapper` overload.
- `gammaphi_gibbs` now handles unnormalized inputs (mole amounts) correctly and always returns a `(g, v)` tuple for empty phases.
