# Estimation basics

To estimate parameters of an equation of state, we need to remember that what we are actually doing is minimizing a loss function $L(\Theta,P) = \sum{L_{i}(\Theta,P_i)}$, where $\Theta$ is a list of parameters (in this particular case, used by an equation of state) and $P$ is one or more datasets.
Once this optimization problem is formulated, we can obtain the optimum list of parameters $\Theta_{opt}$ by using any optimizer.

So, in the particular case of equations of state with Clapeyron.jl (or any other equation of state package), solving an estimation problem requires:

- $\Theta$: a way to switch representations, between a list of arbitrary parameters and an equation of state containing those parameters, in function form, we need `Θ(model)` and `model(Θ)`. In practical terms, $\Theta$ is treated as as vector of numbers.
- $P_i,l_i$: a way to get the loss of a given dataset. In particular, if instead of using the parameters in the loss function ($L_{i}(\Theta,P_i)$), we directly use the equation of state ($L_{i}(EoS(\Theta),P_i)$) we can separate the representantions of datasets and parameters.
- Misc: initial guesses, lower and upper bounds over the parameters.

## the `EstimationModel` struct

This object provides the necessary infraestructure to transform representations, between an `EoSModel` and
a list of parameters. We can build an `EstimationModel` object with an existing equation of state, along with additional data:

```julia
model = PCSAFT("methane")
toestimate = [
    Dict(
        :param => :epsilon,
        :lower => 100.,
        :upper => 300.,
        :guess => 250.
    ),
    Dict(
        :param => :sigma,
        :factor => 1e-10,
        :lower => 3.2,
        :upper => 4.0,
        :guess => 3.7
    )
    ,
    Dict(
        :param => :segment,
        :lower => 0.9,
        :upper => 1.1,
        :guess => 1.
    )
]

est_model = EstimationModel(model,toestimate)
```

With the `EstimationModel` created, we can now interact with the stored equation of state in the following ways:

```julia
import Clapeyron.EstimationUtils
#if you use `using Clapeyron.EstimationUtils` instead, the functions will be loaded into the project instead

#get a vector of parameters from the model
theta = EstimationUtils.parameter_vector(est_model)

#get initial guess
theta0 = EstimationUtils.initial_guess(est_model)

#modify vector of parameters
theta1 = 2 .* theta

#store vector of parameter in EoSModel
EstimationUtils.set_parameter_vector!(est_model,theta1)

#or, we can do it using julia broaccasting:
est_model .= theta1

@assert theta1 == EstimationUtils.parameter_vector(est_model)

lb,ub = EstimationUtils.lower_bounds(est_model),EstimationUtils.upper_bounds(est_model)

#EstimationModel implements the indexing interface for easy access to the parameters
theta2 = 0.7 .* theta
est_model[1] = 0.7*theta[1]
t2 = est_model[2]
```

The information in `toestimate` indicates how our resulting `EstimationModel` will be used, how many parameters does it allow, bounds, initial guesses, and how do the input parameters interact with the stored equation of state.

### Indices

When constructing an `EstimationModel` a mapping is created between the EoS fields and the indices of a parameter vector.
For most models, the mapping is 1-to-1, that is, the index `2` will map to to the model's only `sigma` parameter value.
We can control how this mapping is created via the `indices` keyword.

On single component models, The `indices` keyword can only be ommited if the EoS parameter has only one value.
On multicomponent models, the specification of indices is necessary and useful; we can exploit parameter symmetries to reduce the parameter space, or omit parameters that only depend on single component data.
The index type depends on the type of parameter:

- parameters of type `Clapeyron.SingleParam` are indexed via integers.
- parameters of type `Clapeyron.PairParam` are indexed in two ways: an integer will get the diagonal value of the pair matrix, whereas a tuple of integers is capable to indexing the entire parameter, including off-diagonal entries. By default, we assume symmetry: any changes of the parameter on the entry `(i,j)` will also be done on the entry `(j,i)`. One can turn off this behaviour (on activity models, for example) by passing the keyword `:symmetric => false`.
- parameters of type `Clapeyron.AssocParam` are indexed via integers, in the same way that `SingleParam`s. one index corresponds to each association pair. Some association pairs could have an associated cross-association value, to modify the cross-association value and the specified value simultaneously, we can use the `cross_assoc => true` keyword.
- Some off-diagonal pair parameters are computed as a combination of one or more other EoS parameters that may or may not be used in the current estimation procedure. For example, in the `SAFTVRMie` EoS both `sigma` and `epsilon` are computed via combining rules that depend on both parameters diagonal entries. In particular, `Clapeyron.jl` has a mechanism to allow specific entries of the pair parameter matrix to be marked as *missing* and be recalculated each time a change is done in other parameters of the same model, or to be marked as *not missing*, and be fixed when any recalculation occurs.If one wants to modify diagonal values of a pair parameter that may participate in any recombination procedure, we can use `recombine => true` keyword to mark all entries in the same row and column as the diagonal index as missing, allowing their recombination.

### scaling factors

`Clapeyron.jl` models normally store physically-based parameters in the SI metric system, other model parameters may be stored in arbitrary scales. While the reasons of using one scale or the other depend on the EoS developer, what could be an adequate scale for one application may be not be ideal for other use cases.
For example, the `sigma` parameter on the `PCSAFT` equation of state (the hard-sphere diameter of the monomer) is stored in meters; just to get a sense of scale used, the hard-sphere diameter of methane (used by the PCSAFT EoS) is around `3.7e-10` meters.
While this specific scale is useful when calculating properties that are in the SI metric system, we may prefer another scale when performing parameter estimation.
In particular, we may want to scale all parameters used in the parameter estimation procedure in a way that they all have similar magnitudes, easing the work of the optimizer procedure.
For this particular reason, an scaling `factor` can be specified: When storing a parameter in the equation of state, that parameter is multiplied by the scaling factor; this scaling is reversed when getting a parameter instead, so one can move freely between the scaling used in the estimation procedure and the scaling used by the equation of state.

### Bounds, guesses

In addition to information on how to access and modify an equation of state, most optimizers will require lower and upper bounds for each parameter. Some optimizers will additionally require an initial guess from which start their optimization procedure. Bounds can be provided by the `lower` and `upper` keywords, while the initial guess is provided by the `guess` keyword. The default bounds are the entire real line (`(-Inf,Inf)`), the initial guess value is simply extracted from the equation of state if not specified.

## the `EstimationData` struct

The `EstimationData` object stores the following properties:

- a set of input data points (and their errors if available)
- a set of output data points (and their errors if available)
- weights for each data point (by default all points have weight equal to `1`)
- a method `f` that takes a determinate amount of inputs and returns a determinate amount of outputs.
- a loss function `L` that takes two terms: the result of `f(input)` and the expected `output`.

The `EstimationData` constructor takes a CSV file or other column-based table (only one file per `EstimationData` object) and it sorts the columns according the following criteria:

- Column names starting with `out_` are considered output values
- Column names that dont follow any of this criteria are considered input calues

```julia

function my_sat_pressure(model,T)
    saturation_pressure(model,T)[1]
end

Tsat = [90.694, 100.69, 110.69, 120.69, 130.69, 140.69, 150.69, 160.69, 170.69, 180.69]
p_sat = [11696.0, 36936.0, 93451.0, 201010.0, 382860.0, 664510.0, 1.073e6, 1.6369e6, 2.3872e6, 3.361e6]
psat_data = (T = Tsat,out_p = p_sat)
my_loss(x,y) = abs2(x - y)
est_data = EstimationData(my_sat_pressure,psat_data,my_loss)
```

We can now use our `EstimatorData` object to evaluate the loss function:

```julia
import Clapeyron.EstimationUtils
model1 = cPR("methane")
model2 = SingleFluid("methane")

loss1 = EstimationUtils.objective_function(est_data,model1)
loss2 = EstimationUtils.objective_function(est_data,model2)
```

###

## the `EstimationProblem` struct

With an `EstimationModel` object containing