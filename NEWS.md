# v0.6.25

## New Features

- Estimation framework: revamped internals, what was previously `Estimation` was renamed `EstimationProblem`, and there are two new structs: `EstimationModel`, a wrapper over a model and `ToEstimate` that allows indexing models as if they are the vector of parameters described by `ToEstimate`; and `EstimationData`, that is a reimplementation of the old `EstimationData` struct that additionally stores the loss and the estimator method, allowing multiple loss functions per estimation problem, along with standalone, model-independent evaluations of the objective function. Current estimation workflows that used the documented API are not affected, apart from performance improvements.
To see a full list of the changes, we recomend checking the new estimation basics tutorial.
- Activity models: a consistent gibbs bulk property model for the liquid phase was defined, based on the incompressible approximation already used on Clapeyron equilibria solvers. The gas phase is still evaluated in the helmholtz framework.
- New flash method: `RRXYFlash`, a X-(T or P) flash with support for activity and composite models.
- Flash: activity and composite models can now get properties via `property(model,result::FlashResult)`
- Flash: `FlashResult` now stores the index of vapour phase. it can be accessed via `identify_phase(result,i)`
- Activity models: support for `Tproperty`/`Pproperty`/`edge_pressure`/`edge_temperature`
- New model: Patel-Teja-Hayen cubic (`PatelTejaHayen`)
- New alpha model: Twu alpha estimation used in `tcPR` and `tcRK` (`tcTwuAlpha`)
- New translation model: `tcPR` and `tcRK` translation estimation (`tcTranslation`)
- Cubic mixing rules: support for more cubics. Before, most EoS + GE mixing rules were only correct for cubics with composition-independent coefficients.
- CSV parsing: options can now be provided via an inline JSON file: `[csvtype = single]` is equivalent to `{"csvtype" : "single"}`
- Flash methods: added `verbose` option to `RRQXFlash`, improved `verbose` display for `RRTPFlash`/`MichelsenTPFlash`
- Various bug fixes and stability improvements.

## Other Changes

- `ClapeyronHANNA` was deprecated. the new version just uses the implementation on `MLThermoProperties` instead.
