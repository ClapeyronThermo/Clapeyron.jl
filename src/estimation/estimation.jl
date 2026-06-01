include("estimation_utils.jl")
include("estimation_data.jl")
include("estimation_model.jl")

"""
    EstimationProblem(est_model::AbstractEstimationModel,data;concrete = false)

Core structure used for parameter optimization. 
It joins estimation models and estimation data to perform parameter optimization.
It can be created from a `EstimationUtils.AbstractEstimationModel` and a list of `EstimationUtils.AbstractEstimationLoss`.
If `concrete` is set to `true`, then the list of data will be converted into a tuple before storing it.
If `concrete` is set to `false`, then the list of data will be stored as an abstract vector, allowing adding data with different losses and methods afer the problem is constructed.

"""
mutable struct EstimationProblem{T<:EoSModel, M <: EstimationUtils.AbstractEstimationModel{T},D}
    model::T #this model is an alias of the model stored inside toestimate.
    initial_model::T #we dont touch this particular model
    toestimate::M
    data::D #abstractly typed for easy update
end
# Mutable for now to make it easy to just replace the model

function EstimationProblem(est_model::EstimationUtils.AbstractEstimationModel,data;concrete = true)
    if concrete
        new_data = tuple(data...)
    else
        new_data = Vector{EstimationUtils.AbstractEstimationLoss}[]
        resize!(new_data,length(data))
        new_data .= data
    end

    model = EstimationUtils.get_model(est_model)
    model2 = deepcopy(model)
    est_model2 = EstimationUtils.set_model(est_model,model2)

    comps = component_list(model)
    norm_comps = normalisestring.(comps)
    for data_i in new_data
        __estimationdata_fix_species!(model,data_i,comps,norm_comps)
    end

    return EstimationProblem(model2,model,est_model2,new_data)
end

__estimationdata_fix_species!(data,comps,norm_comps) = nothing

function __estimationdata_fix_species!(data::EstimationData,comps,norm_comps)
    #additional checks
    species = data.species
    species[1] == "all" && length(species) == 1 && return nothing
    for (j,species_j) in enumerate(species)
        if species_j ∉ comps
            norm_species_j = normalisestring(species_j)
            i = findfirt(isequal(norm_species_j),norm_comps)
            if isnothing(i)
                throw(error("EstimationData error: species $species_j not found in input model, please check the input `species` field in the CSV."))
            else
                #modify species to make it match with the input model
                species[j] = comps[i]
            end
        end
    end
end

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
function Estimation end

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

function Estimation(est_model::EstimationUtils.AbstractEstimationModel,data)
    return EstimationProblem(est_model,toestimate)
end

function _Estimation(model::EoSModel, toestimate::Vector{Dict{Symbol,Any}}, filepaths, objective_form, ignorefield) 
    est_model = EstimationModel(model,toestimate,ignorefield)
    estimation = EstimationProblem(est_model, estimation_data_from_csvs(filepaths, objective_form),concrete = false)
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

function reload_data(estimation::EstimationProblem)
    for (i,est_data) in enumerate(estimation.data)
        source = est_data.source
        _loss = est_data_options.loss == Symbol() ? est_data.loss : nothing
        _method = est_data_options.method == Symbol() ? est_data.method : nothing
        _output_weights = est_data_options.output_weights == [1.0] ? est_data.output_weights : 1.0
        push!(method,EstimationData)
        estimation.data[i] = EstimationData((_output_weights,source),_method,_loss)
    end
    return estimation
end

#update via Symbols
function update_estimation!(estimation::EstimationProblem, params::Vector{Symbol},  values) 
    indices = EstimationUtils.symbol_indices(estimation)
    @assert length(indices) == length(values)
    Θ₀ = EstimationUtils.get_eos_parameters(estimation)
    Θ = @view Θ₀[indices]
    Θ .= values
    set_eos_parameters!(estimation,Θ)
    return estimation
end

#update model
function update_estimation!(estimation::EstimationProblem, model::EoSModel)
    estimation.model = model
    return estimation
end

#update data
function update_estimation!(estimation::EstimationProblem, filepaths::Union{AbstractString,AbstractVector{<:AbstractString}})
    empty!(estimation.data)
    new_data = estimation_data_from_csvs(filepaths)
    resize!(estimation_data,length(new_data))
    estimation.data .= new_data
    return estimation
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
        eval_funcion_i = estimation_data[i]
        Fi = objective_function(eval_funcion_i,Θmodel)
        if !isfinite(Fi)
            F = 1e100*oneunit(F)
            break
        end
        F += Fi
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
EstimationUtils.get_eos_parameters(model::EstimationProblem) = EstimationUtils.get_eos_parameters(model.toestimate)

export Estimation, EstimationModel, EstimationData, ToEstimate
