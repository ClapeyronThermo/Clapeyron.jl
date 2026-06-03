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
theta = EstimationUtils.get_eos_parameters(est_model)

#get initial guess
theta0 = EstimationUtils.initial_guess(est_model)

#modify vector of parameters
theta1 = 2 .* theta

#store vector of parameter in EoSModel
EstimationUtils.set_eos_parameters!(est_model,theta1)

#or, we can do it using julia broadcasting:
est_model .= theta1

@assert theta1 == EstimationUtils.get_eos_parameters(est_model)

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
- Some off-diagonal pair parameters are computed as a combination of one or more other EoS parameters that may or may not be used in the current estimation procedure. For example, in the `SAFTVRMie` EoS both `sigma` and `epsilon` are computed via combining rules that depend on both parameters diagonal entries. In particular, `Clapeyron.jl` has a mechanism to allow specific entries of the pair parameter matrix to be marked as *missing* and be recalculated each time a change is done in other parameters of the same model, or to be marked as *not missing*, and be fixed when any recalculation occurs. If one wants to modify diagonal values of a pair parameter that may participate in any recombination procedure, we can use `recombine => true` keyword to mark all entries in the same row and column as the diagonal index as missing, allowing their recombination.

### Symbol indices

`EstimationModel` also supports indexing via `Symbol`, with the slight caveat that `est_model[sym]` returns a vector of values instead of a single value:

```julia
sigma = est_model[:sigma]
@assert sigma isa Vector #true
est_model[:sigma] = 3 #works
est_model[[:sigma,:epsilon]] = [3,300] #works

#limitation: this does not work, it will not set the indices.
est_model[[:sigma,:epsilon]] .= [3,300] 
```

### Multiple-value indices

`EstimationModel` also supports multiple indices per dictionary:

```julia
model = PCSAFT(["water","ethanol"])

toestimate = [
    Dict(
        :param => :epsilon,
        :indices => [1,2], #default to epsilon[1,1] and epsilon[2,2]
        :lower => 100., #defaults to [100,100]
        :upper => [300.,400],
        :guess => 250. #defaults to [250, 250]
    ),
    Dict(
        :param => :sigma,
        :indices => [(1,1),(1,2)], #pairs of values are also accepted
        :factor => 1e-10,
        :lower => 3.2,
        :upper => 4.0,
        :guess => 3.7
    )
]

est_model_multiple = EstimationModel(model,toestimate)
EstimationUtils.parameter_length(est_model_multiple) #4
```

### Special multiple-value indices

There are also some special indices that can be used by passing the corresponding symbol to the `indices` key:

- `:pures` — includes only the pure-component (diagonal) entries. For `PairParam` it selects `(i,i)` for each component `i`; for `AssocParam` it selects all association pairs where the two site indices belong to the same component.
- `:unlike` — includes only the cross (off-diagonal) entries. For `PairParam` these are the entries `(i,j)` with `i ≠ j`; for `AssocParam` these are association pairs that span two different components. On `SingleParam` this special index will fail.
- `:all` — includes every entry, diagonal and off-diagonal.

!!! note

    On earlier package versions, if no `indices` key is specified and the model has more than one component, the default is set to the first index of the respective parameter. This behaviour could change on future versions of `Clapeyron.jl`

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

- Column names starting with `out_` are considered output values.
- Column names ending in `_error_abs`,`error_rel`,`error_std` are considered error values. Those are not utilized in the parameter estimation procedure, but they are useful for calculating statistical properties of the fitted parameters later.
- Column names that dont follow any of this criteria are considered input values.

```julia
#Experimental Melting and Vapor Pressures of Methane, J. Chem. Thermodyn., 1972, 4, 1, 127-133, https://doi.org/10.1016/S0021-9614(72)80016-8 .
A,B,C = (3.9895,443.028,-0.49)
antoine_psat(T) = exp10(A - (B/(T + C)))*1e5
Tsat = collect(range(90.99,185.99,step = 5))
psat = antoine_psat.(Tsat)

#compiled data
psat_data = (T = Tsat, out_p = psat)

#loss function
my_loss(y_calc,y_exp) = abs2(y_calc - y_exp)

#function to get evaluated properties:
my_sat_pressure(model,T) = saturation_pressure(model,T)[1]

est_data = EstimationData(my_sat_pressure,psat_data,my_loss)
```

We can now use our `EstimationData` object to evaluate the loss function:

```julia
import Clapeyron.EstimationUtils
model1 = cPR("methane")
model2 = SingleFluid("methane")

loss1 = EstimationUtils.objective_function(est_data,model1)
loss2 = EstimationUtils.objective_function(est_data,model2)
```

### Multiple function returns

Most phase equilibria procedures don't return an unique value. Instead, they return all information necessary to reconstruct the phases.
For example, `Clapeyron.saturation_pressure(model,T)` returns not only the saturation pressure, but also the liquid and vapour volumes at `(Psat,Tsat)` conditions.
While performing a parameter estimation procedure, it is common to find that two or more loss functions depend on only one phase equilibria procedure.
For example,the common parameter estimation procedure for the PC-SAFT EoS for pure, non-associating molecules is to minimize the loss with respect to the saturated liquid density and saturation pressure, but both properties are the same thermodynamic state: pure saturated liquid at (Psat,Tsat), and can be obtained simultaneusly.

The `EstimationData` struct accounts for this automatically. One just needs to provide multiple output values:

```julia

#same temps as before
A,B,C = (3.9895,443.028,-0.49)
antoine_psat(T) = exp10(A - (B/(T + C)))*1e5
Tsat = collect(range(90.99,185.99,step = 5))
psat = antoine_psat.(Tsat)

#we are gonna extract liquid density values from the reference equation of state for methane
truth_model = SingleFluid("methane")
Vlsat = first.(Clapeyron.x0_sat_pure.(truth_model,Tsat)) #we are gonna use the ancillary for the liquid volume values.

#fitting density instead of volume:
rholsat = 1 ./ Vlsat

#compiled data with both properties
psat_and_rholsat_data = (T = Tsat, out_p = psat, out_rholsat = rholsat)

#our method will be different:
function psat_and_rhosat(model,T)
    p,vl,_ = saturation_pressure(model,T)
    return p,1/vl
end

#the loss will be applied to (p_sat - p_exp) and (rhol - rhol_exp)
my_loss(y_calc,y_exp) = abs2(y_calc - y_exp)

est_data = EstimationData(psat_and_rhosat,psat_and_rholsat_data,my_loss)
```

On functions with multiple return values, we just sum the losses of each output value with their respective expected data points. While this is ok for some simple properties, we may want to add a weight to each output value; some properties may have more error, or they may be just estimations from correlations.
`EstimationData` supports adding weighting values via the `output_weights` keyword argument. Using the same example as before:

```julia
#now, the loss of each data point will be equal to 4*loss(psat - psat_exp) + 1.1*loss(rhol - rhol_exp)
est_data = EstimationData(psat_and_rhosat,psat_and_rholsat_data,my_loss,output_weights = (4.0,1.1))
```

### Multiple function inputs

Conversely, in the same way that could be multiple properties per data point, There are properties that require multiple inputs.
Most single-component data away from saturation depends on the pressure *and* the temperature, multiple component models also depend on the composition of the mixture.
`EstimationData` handles this situation similar to how its handled with multiple outputs, we just need to pass more input vectors:

```julia
#fitting density of a gas:

Tx = 400:5:500
px = 1e4:1e4:1e5

#all combinations of (px,Tx)
ptx = vec(collect(Iterators.product(px,Tx)))
p = first.(ptx)
T = last.(ptx)

#IAPWS-95 as reference data
truth_model = IAPWS95()
rhov = mass_density.(truth_model,p,T,phase = :v)

#p first then T, matching the argument order of my_den below
rhov_data = (p = p, T = T, out_rhov = rhov)

#multiple input parameters, in the same order as the table
my_den(model,p,T) = mass_density(model,p,T,phase = :v)

my_loss(y_calc,y_exp) = abs2(y_calc - y_exp)

est_data = EstimationData(my_den,rhov_data,my_loss)
```

### Reading from a CSV

instead of defining everything in code, it is useful to store, along with the dataset, the expected loss functions, evaluation methods, among other properties. `EstimationData` is capable of reading a CSV with the following format:

```csv
Clapeyron Estimator #main header 
[method = my_method,loss = my_loss,normalize = true,output_weights = 1.0 2.0 3.0,data_weights_col = dw,species = "carbon dioxide" "methane",]
x,y,z,out_1,out_2,out_3
```

When reading a CSV, `EstimationData` parses the options line (the second line, enclosed in `[...]`) to extract the method, loss, weights, and other settings. Any option that is also provided as a keyword argument to the constructor takes priority over the value found in the CSV. If neither the CSV nor the keyword argument specifies a loss, the default mean squared relative error ($L(x,y) = (\frac{x - y}{y})^2$) is used. If no method is found, an error is raised.

The `method` and `loss` keys are looked up by name in the `Main` module, so they must be defined there before constructing the `EstimationData`. Species names are parsed as a space-separated list of quoted strings; if not provided, all components of the model are used.

```julia
#given a CSV file "psat_methane.csv" with a matching method and loss already defined in Main:
est_data = EstimationData("psat_methane.csv")

#keyword arguments override the CSV options:
est_data = EstimationData("psat_methane.csv", my_sat_pressure, my_loss; normalize = false)
```

## the `EstimationProblem` struct

With an `EstimationModel` object containing the parameter vector and an `EstimationData` object (or a collection of them) containing the experimental data, we have everything needed to define a full estimation problem. The `EstimationProblem` struct bundles both together and exposes a single `objective_function` that the optimizer can call:

```julia
#using the est_model and est_data defined earlier
prob = EstimationProblem(est_model, est_data)
```

`EstimationProblem` also accepts a vector of `EstimationData` objects when fitting against multiple datasets simultaneously. The total loss is the sum of the individual losses:

```julia
est_data_psat  = EstimationData(my_sat_pressure, psat_data, my_loss)
est_data_rhol  = EstimationData(my_den, rhov_data, my_loss)

prob = EstimationProblem(est_model, [est_data_psat, est_data_rhol])
```

To evaluate the objective function at a given parameter vector `Θ`:

```julia
import Clapeyron.EstimationUtils

Θ0 = EstimationUtils.initial_guess(prob)
lb = EstimationUtils.lower_bounds(prob)
ub = EstimationUtils.upper_bounds(prob)

loss = EstimationUtils.objective_function(prob, Θ0)
```

Internally, `objective_function(prob, Θ)` calls `set_eos_parameters!(est_model, Θ)` to update the wrapped EoS model, then sums `objective_function(data_i, est_model.model)` over all `EstimationData` objects.

### Running the optimization

Once an `EstimationProblem` is assembled, We can pass it to any optimizer to solve our estimation problem. A common choice is a gradient-free box-constrained method, such as those available in `Optim.jl`:

```julia
using Optim

Θ0 = EstimationUtils.initial_guess(prob)
lb = EstimationUtils.lower_bounds(prob)
ub = EstimationUtils.upper_bounds(prob)
f = EstimationUtils.objective_function(prob) #equivalent to Θ -> EstimationUtils.objective_function(prob, Θ)
result = optimize(
    f,
    lb, ub, Θ0,
    Fminbox(NelderMead())
)

Θ_opt = Optim.minimizer(result)
```

On the choice of optimizers and julia packages, We have been obtaining good fits with the optimizer methods found in `Metaheuristics.jl`. There is also explicit support via package extensions, to pass the `EstimationProblem` directly to the optimizer:

```julia
using Metaheuristics
Θ_opt, fitted_model = Metaheuristics.optimize(prob,ECA(;options=Options(iterations=100)))
```

### Retrieving the solution

After the optimization is complete, you can apply the optimal parameters back to the model and inspect the result:

```julia
#write the optimal parameters into the EoS model
EstimationUtils.set_eos_parameters!(est_model, Θ_opt)

#the underlying EoS model now uses the fitted parameters
fitted_model = EstimationUtils.get_model(est_model)

#compute any property with the fitted model as usual
p_fit, vl_fit, vv_fit = saturation_pressure(fitted_model, 150.0)
```

The fitted model is a standard `EoSModel` and can be used everywhere a regular Clapeyron model is accepted, including property calculations, phase diagrams, and further refinement.
