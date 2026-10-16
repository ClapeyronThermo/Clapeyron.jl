```@meta
CurrentModule = Clapeyron
```

## Contents

```@contents
Pages = ["estimation.md"]
```

## Index

```@index
Pages = ["estimation.md"]
```

## Estimation core API

The estimation API is stored inside the `EstimationUtils` submodule.

```@docs
Clapeyron.EstimationUtils
```

### `AbstractEstimationModel` API

```@docs
Clapeyron.EstimationUtils.AbstractEstimationModel
Clapeyron.EstimationUtils.get_eos_parameters
Clapeyron.EstimationUtils.set_eos_parameters!
Clapeyron.EstimationUtils.get_model
Clapeyron.EstimationUtils.set_model
Clapeyron.EstimationUtils.parameter_length
Clapeyron.EstimationUtils.symbol_indices
Clapeyron.EstimationUtils.lower_bounds
Clapeyron.EstimationUtils.upper_bounds
Clapeyron.EstimationUtils.initial_guess
```

### `AbstractEstimationLoss` API

```@docs
Clapeyron.EstimationUtils.AbstractEstimationLoss
Clapeyron.EstimationUtils.objective_function
```

## Estimation implementation

We implement the Estimation API defined in the `EstimationUtils` module via two structures:

- `EstimationModel` implements the `AbstractEstimationModel` API
- `EstimationData` implements the `AbstractEstimationLoss` API

along with those structures, we implement the `EstimationProblem` struct, that is a wrapper over an `AbstractEstimationModel` and a list of `AbstractEstimationLoss`, with support for some global optimizers.

```@docs
Clapeyron.EstimationData
Clapeyron.EstimationModel
Clapeyron.EstimationProblem
Clapeyron.Estimation
```

## Estimation utilities

```@docs
Clapeyron.reload_data!
```
