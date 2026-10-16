# Custom estimation types

The built-in [`EstimationModel`](@ref) and [`EstimationData`](@ref) structs cover the most common parameter estimation workflows in Clapeyron.jl, but there are situations where a more direct implementation is preferable: specialised parameter layouts, non-standard property calculations, or simply avoiding the overhead of the general machinery. This guide shows how to implement the two extension points of the framework from scratch.

- **`AbstractEstimationModel{M}`** — controls how a flat parameter vector `Θ` maps to and from an EoS model. The example fits binary interaction parameters `kij` via `Clapeyron.get_k` / `Clapeyron.set_k!`.
- **`AbstractEstimationLoss`** — computes the scalar objective function from experimental data. The example fits saturation pressure data directly, without going through the CSV machinery of `EstimationData`.

## Custom `AbstractEstimationModel`: fitting `kij` parameters of binary cubic models

In the specific case of cubic models, binary interaction parameters (BIPs) are indirectly stored in the off-diagonal values of `a` and `b`. how they are stored depends on the specific mixing rule.

The required interface for `AbstractEstimationModel{M}` is:

| Function | Purpose |
|---|---|
| `EstimationUtils.get_eos_parameters(est)` | Return the flat vector `Θ` extracted from the wrapped model |
| `EstimationUtils.set_eos_parameters!(est, Θ)` | Write `Θ` back into the wrapped model in-place |
| `EstimationUtils.get_model(est)` | Return the wrapped EoS model |
| `EstimationUtils.set_model(est, new_model)` | Return a new wrapper around `new_model` |

Optionally, we can implement `lower_bounds`, `upper_bounds`, and `initial_guess` so that the struct carries all the information an optimiser needs without extra arguments.

### Struct definition and API implementation

```julia
import Clapeyron
import Clapeyron.EstimationUtils
using StaticArrays

struct KijEstimationModel{M<:Clapeyron.CubicModel} <: EstimationUtils.AbstractEstimationModel{M}
    model::M
end

#AbstractEstimationModel mandatory interface
function EstimationUtils.get_eos_parameters(est_model::KijEstimationModel)
    [Clapeyron.get_k(est_model)[1,2]]
end

function EstimationUtils.set_eos_parameters!(est_model::KijEstimationModel,k12)
    k = @SVector [k12;zero(k12);;zero(k12);k12]
    Clapeyron.set_k!(k)
end

#we dont need to define EstimationUtils.set_model or EstimationUtils.get_model because we sotre our model in the `model` field.

#helper functions

#kij = 0 by default
EstimationUtils.initial_guess(est_model::KijEstimationModel) = [zero(eltype(est_model.model))]

#=at kij = 1, we have (1 - kij) = 0
but this is not the theoretical upper bound.
#for a positive a value, the upper bound is dependent on the composition:

k12 < 1 + (x1^2 * a11 + x2^2 * x22)/(2*x1*x2*sqrt(a11*a22))
=#
EstimationUtils.upper_bounds(est_model::KijEstimationModel) = [oneunit(eltype(est_model.model))]


#just an arbitrary lower value
EstimationUtils.loewer_bounds(est_model::KijEstimationModel) = [-oneunit(eltype(est_model.model))]
```

We can now set/get kij for binary models via the `AbstractEstimationModel` API:

```julia
# Two-component PR model for CO₂ + methane
model = PR(["carbon dioxide", "methane"])

# Wrap it; no explicit guess → use the current kij values (typically 0.0)
est_model = KijEstimationModel(model)

# Read the current kij 
Θ0 = EstimationUtils.initial_guess(est_model)   # [0.0]

# Write a new value
EstimationUtils.set_eos_parameters!(est_model, [0.1])

# The underlying model is updated in-place
@assert Clapeyron.get_k(est_model.model)[1, 2] ≈ 0.1
```

## Custom `AbstractEstimationLoss`: saturation pressure fitting with reference model

Instead of building an `EstimationData` from a CSV, we can encode the experimental data and the property calculation directly inside a custom `AbstractEstimationLoss` subtype. This is useful when the data is generated programmatically (e.g. from a reference EoS) or when the loss calculation needs logic that does not fit neatly into a single method function.

The required interface for `AbstractEstimationLoss` is a single function:

| Function | Purpose |
|---|---|
| `EstimationUtils.objective_function(loss, model)` | Evaluate the scalar loss for `model` against the stored data |

```julia
struct PsatWithReference{REF} <: EstimationUtils.AbstractEstimationLoss
    reference_model::REF
    T_exp::Vector{Float64}
    P_cache::Dict{Float64,Float64}
end

#convenience constructor
PsatWithReference(name::String,T_exp) = PsatWithReference(SingleFluid("methane"),T_exp,Dict{Float64,Float64}())

function EstimationUtils.objective_function(loss::PsatWithReference, model)
    T_exp = loss.T_exp
    P_cache = loss.P_cache
    reference_model =  loss.reference_model
    # If model is a mixture, isolate the component we care about.

    reference_model_name = Clapeyron.component_list(reference_model)[1]
    model_names = Clapeyron.component_list(model)
    
    if reference_model_name in model_names
        model_r = if length(model) == 1
                model
            else
                ix = findfirst(isequal(reference_model_name),model_names)
                Clapeyron.each_split_model(model,[i])
            end
    else
        error("species $reference_model_name not found in model")
    end

    F   = zero(eltype(model))
    valid = 0

    for T in T_exp
        p_ref = get!(P_cache,T_exp) do
            saturation_pressure(reference_model,T_exp)[1]
        end

        if isfinite(p_calc) && isfinite(p_ref)
            F     += abs2((p_calc - p_ref) / p_ref)
            valid += 1
        end
    end

    return valid > 0 ? F / valid : Inf*one(typeof(F))
end
```

The function returns the mean squared relative error $\frac{1}{N}\sum_{i=1}^{N}\left(\frac{p_\mathrm{calc} - p_\mathrm{exp}}{p_\mathrm{exp}}\right)^2$, normalised by the number of valid (finite, positive) data points, and `Inf` if none are valid.

We can now use it

```julia
# Generate reference saturation pressures from a high-accuracy model.
ref_model = SingleFluid("methane")
T_exp = collect(range(100.0, 180.0, step = 5.0))

psat_loss = PsatWithReference(ref_model,T_exp,Dict{Float64,Float64}())

# Evaluate against a candidate model.
candidate = PCSAFT("methane")
loss_val = EstimationUtils.objective_function(psat_loss, candidate)
```

---

## Putting it together: fitting kij from binary VLE data

With both custom types in hand, we can now assemble a full estimation problem. For a `KijEstimationModel` we use `EstimationProblem` directly; `PsatLoss` is used here to show that any `AbstractEstimationLoss` can be composed into the same workflow.

```julia
#model to fit
model    = PR(["carbon dioxide", "methane"])
est_model = KijEstimationModel(model)

#functions used in the CSVs
function dew_point_fn(model,y1,p) 
    T,_,_,x = dew_temperature(model,p,[y1,1-y1])
    return p,x[1]
end

function bubble_point_fn(model,x1,p) 
    T,_,_,y = dew_temperature(model,p,[x1,1-x1])
    return T,y[1]
end

data_sources = ["bubble_point.csv","dew_point.csv"]
#read CSV via utility fn
data = Clapeyron.estimationdata_from_csvs(data_sources)

#Estimation problem
prob = EstimationProblem(est_model, data)
f = EstimationUtils.objective_function(prob)
#using Metaheuristics directly
Θ_opt, _  = Metaheuristics.optimize(prob, ECA());

# Write the optimal kij back into the model.
EstimationUtils.set_eos_parameters!(est_model, Θ_opt)
fitted_model = EstimationUtils.get_model(est_model)

println("Optimal kij: ", Θ_opt[1])
```

### Composing multiple loss terms

Because `EstimationProblem` accepts a vector of `AbstractEstimationLoss` objects, we can combine our custom `PsatWithReference` with a standard `EstimationData` (or another custom loss) and the total objective is their sum:

```julia
# Saturation pressure loss from the custom type.
psat_loss = PsatWithReference("carbon dioxide", T_exp)

# Liquid density loss from the standard EstimationData machinery.
rhol_method(model, T) = 1.0 / Clapeyron.saturation_pressure(model, T)[2]
rhol_data  = (T = T_exp, out_rhol = rhol_exp)
rhol_loss  = EstimationData(rhol_method, rhol_data, __mse)

# Both losses are minimised simultaneously.
prob = EstimationProblem(est_model, [psat_loss, rhol_loss])
```
