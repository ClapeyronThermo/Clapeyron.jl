include("estimation_utils.jl")
include("estimation_data.jl")
include("estimation_model.jl")

"""
    est,obj,x0,lb,ub = Estimation(model::EoSModel,toestimate::Dict,filepaths;ignorefield = Vector{String},objective_form = mse(pred,exp) = ((pred-exp)/exp)^2)
    est,obj,x0,lb,ub = Estimation(est_model::EstimationModel,to_estimate::Union{Vector{EstimationData},NTuple{N,EstimationData}) where N
    est,obj,x0,lb,ub = Estimation(est_model::EstimationModel,to_estimate::EstimationData)

## Input parameters:
- `model`: The initial model containing the species we wish to parameterise
- `toestimate`: The dictionary of parameters being fitted, or already instantiated `EstimationData` objects.
- `filepaths` or `filepaths_weights`: The location of the data files used to fit. Can also contain the weights of each dataset
- `ignorefield`: Specify which EoSModel fields to ignore in the main model
- `objective_form`: Specify the functional form of the objective function in the form `objective_form(pred,exp)`

## Output:

- `est`: an `EstimationProblem` object which contains the following fields:
    - `model`: The model whose parameters will be varied
    - `initial_model`: The initial model before parameterisation
    - `toestimate`: `EstimateModel` struct, which contains information on how to transform between an specified model and a vector of parameters
    - `data`: a collection of `EstimationData` objects. each one contributing to the objective function `obj`
- `obj`: The objective function which is used to fit the parameters. It can be also be created via `Clapeyron.EstimationUtils.objective_function(est)`
- `x0`: Initial guesses for the parameters. They can also be accessed via `Clapeyron.EstimationUtils.initial_guess(est)`
- `ub`: Upper bounds for the parameters. They can also be accessed via `Clapeyron.EstimationUtils.upper_bounds(est)`
- `lb`: Lower bounds for the parameters. They can also be accessed via `Clapeyron.EstimationUtils.upper_bounds(est)`

## Description
Produces the estimator and other useful objects used within parameter estimation
"""
mutable struct EstimationProblem{T<:EoSModel,DD}
    model::T
    initial_model::T
    toestimate::EstimationModel{T}
    data::DD
end
# Mutable for now to make it easy to just replace the model

__mse(pred,exp) = ((pred-exp)/exp)^2

function Base.show(io::IO, mime::MIME"text/plain", estimation::EstimationProblem)
    print(io, typeof(estimation))
    println(io, " with data for:")
    show_pairs(io,(d.method for d in estimation.data),prekey = "  :",quote_string = false)
    println(io, "\n to estimate:")
    function val_print(io,val)
        if val !== nothing
            print(io, " with indices => ")
            print('[')
            show_pairs(io,val,nothing,quote_string = false,pair_separator = ',')
            print(']')
        end
    end
    show_pairs(io,estimation.toestimate.params,estimation.toestimate.indices," ",val_print,quote_string = false,prekey = "  :")
end

function Base.show(io::IO, estimation::EstimationProblem)
    print(io, typeof(estimation))
end

function Estimation(model::EoSModel, toestimate::Vector{Dict{Symbol,Any}}, filepaths::Union{Array{String},Array{Tuple{Float64, String}}},objective_form::Base.Callable,ignorefield::Union{Symbol,Vector{Symbol}})
    return _Estimation(model, toestimate, filepaths, objective_form, ignorefield)
end

function Estimation(model::EoSModel, toestimate::Vector{Dict{Symbol,Any}}, filepaths::Union{Array{String},Array{Tuple{Float64, String}}}, ignorefield::Union{Symbol,Vector{Symbol}})
    return _Estimation(model, toestimate, filepaths, __mse, ignorefield)
end

function Estimation(model::EoSModel, toestimate::Vector{Dict{Symbol,Any}}, filepaths::Union{Array{String},Array{Tuple{Float64, String}}},objective_form::Base.Callable)
    return _Estimation(model, toestimate, filepaths, objective_form, Symbol[])
end

function Estimation(model::EoSModel, toestimate::Vector{Dict{Symbol,Any}}, filepaths::Union{Array{String},Array{Tuple{Float64, String}}}; objective_form = __mse, ignorefield = Symbol[])
    return _Estimation(model, toestimate, filepaths, objective_form, ignorefield)
end

function Estimation(est_model::EstimationModel,toestimate::Union{EstimationData,AbstractVector,Tuple})
    if toestimate isa AbstractVector || to_estimate isa Tuple
        @assert all(x -> x isa EstimationData,toestimate)
    end
    concrete_toestimate = tuple(toestimate...)
    model = est_model.model
    est = EstimationProblem(model,deepcopy(model),est_model,concrete_toestimate)
end

function _Estimation(model::EoSModel, toestimate::Vector{Dict{Symbol,Any}}, filepaths, objective_form, ignorefield) 
    est_model = EstimationModel(model,toestimate,ignorefield)
    estimation = EstimationProblem(model, deepcopy(model), est_model, EstimationData(filepaths), objective_form)
    objective = EstimationUtils.objective_function(estimation)
    x0 = EstimationUtils.initial_guess(estimation)
    upper = EstimationUtils.upper_bounds(estimation)
    lower = EstimationUtils.lower_bounds(estimation)
    return estimation, objective, x0, upper, lower
end

export return_model

"""
    return_model
    return_model(estimation::Estimation,model,params)
    return_model(model::EstimationModel,params)


## Input parameters:
- `estimation`: The estimator object
- `model`: The model whose parameters we are varying
- `params`: The new parameters which we want to change
## Output:
- `model`: The new model with the updated parameters
## Description
Based on the parameters provided and the estimator, a new model is produced from the input.
"""
function return_model(estimation::EstimationProblem,model::EoSModel,values)
    T = Base.promote_eltype(model,values)
    return_model!(estimation,promote_model(T,model),values)
end

function return_model!(
    estimation::EstimationProblem,
    model::EoSModel,
    values)
    return set_eos_parameters!(model,estimation,values)
end

@deprecate return_model! set_eos_parameters!

function reload_data(estimation::EstimationProblem)
    estimationdata = EstimationData(estimation.filepaths)
    empty!(estimation.data)
    for i in 1:length(estimation.filepaths)
        push!(estimation.data, estimationdata[i])
    end
end

export update_estimation!

function update_estimation!(
        estimation::EstimationProblem,
        params::Vector{Symbol},
        values::Vector{Any})

    model = estimation.model
    for (i, param) in enumerate(params)
        current_param = getfield(model.params, param)
        if typeof(current_param) <: SingleParameter
            for (j, value) in enumerate(values[i])
                current_param.values[j] = value
            end
        end
        if typeof(current_param) <: PairParam
            for j = 1:length(model.components)
                for k = 1:length(model.components)
                    current_param.values[j,k] = values[i][j,k]
                end
            end
        end
        if typeof(current_param) <: AssocParam
            if typeof(values[i]) <: Compressed4DMatrix
                for (j, value) in enumerate(values[i].values)
                    current_param.values.values[j] = value
                end
            else
                for (j, value) in enumerate(values[i])
                    current_param.values.values[j] = value
                end
            end
        end
    end
end

function update_estimation!(estimation::EstimationProblem, model::EoSModel)
    estimation.model = model
end

"""
    objective_function
    objective_function(estimation::EstimationProblem,params)
    objective_function(estimation::EstimationData,model::EoSModel)
    objective_function(estimation::EstimationProblem)

## Input parameters:
- `estimation`: The estimator object
- `params`: The new parameters which we want to evaluate the objective function for

## Output: The related loss function given the data and parameters provided

## Description
The objective function used within parameter estimation.
"""
function EstimationUtils.objective_function(estimation::EstimationProblem{M,DD},Θ) where {M,DD}
    Θmodel = return_model(estimation, estimation.model, Θ)
    F = zero(Base.promote_eltype(Θmodel))
    data = estimation.data
    for i ∈ 1:length(data)
        F += objective_function(data[i],Θmodel)
    end
    return F
end

function EstimationUtils.objective_function(estimation::EstimationProblem{M,DD},Θ) where {M,DD<:EstimationData}
    Θmodel = return_model(estimation, estimation.model, Θ)
    return objective_function(estimation.data,Θmodel)
end

function EstimationUtils.objective_function(estimation::EstimationProblem)
    return Base.Fix1(EstimationUtils.objective_function,estimation)
end

EstimationUtils.lower_bounds(model::EstimationProblem) = EstimationUtils.lower_bounds(model.toestimate)
EstimationUtils.upper_bounds(model::EstimationProblem) = EstimationUtils.upper_bounds(model.toestimate)
EstimationUtils.initial_guess(model::EstimationProblem) = EstimationUtils.initial_guess(model.toestimate)
EstimationUtils.parameter_vector(model::EstimationProblem) = EstimationUtils.parameter_vector(model.toestimate)

export Estimation, EstimationModel, EstimationData, ToEstimate
