# Estimation basics

To estimate parameters of an equation of state, we need to remember that what we are actually doing is minimizing a loss function $L(\Theta,P) = \sum{L_{i}(\Theta,P_i)}$, where $\Theta$ is a list of parameters (in this particular case, used by an equation of state) and $P$ is one or more sets of data.
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

est_model = EstimationModel(model,data)
```

Some observations:

- Each `Dict` specifies a parameter type, along with lower and upper bounds, and an initial guess. if bounds are not specified, then they are set to `-Inf/Inf`. if the initial guess is not specified, then it is extracted from the model used as an input.
- The dictionary with the `sigma` parameter has an additional `factor` value. This factor is used to convert between the input numeric value (in angstroms) and the numeric value used by the equation of state (meters)
- If a guess is not specified, then it will be taken from the input equation of state model.

Now, we can interact with the equation of state in the following ways:

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

#if we only have one parameter per name, we can use index access via getindex,setindex!

theta2 = 0.7 .* theta
est_model[1] = 0.7*theta[1]
t2 = est_model[2]
t3 = est_model[:segment] #using a symbol instead of a number also works
```

## the `EstimationData` struct

The `EstimationData` object stores the following properties:

- a set of input data points (and their errors if available)
- a set of output data points (and their errors if available)
- weights for each data point (by default all points have weight equal to `1`)
- a method `f` that takes a determinate amount of inputs and returns a determinate amount of outputs.
- a loss function `L` that takes two terms: the result of `f(input)` and the expected `output`.

The `EstimationData` constructor takes a CSV file or other column-based table (only one file per `EstimationData` object) and it sorts the columns according the following criteria:

- Column names starting with `out_` are considered output values
- 
- Column names that dont follow any of this criteria are considered input calues
 

```julia

data = 



```

## the `EstimationProblem` struct
