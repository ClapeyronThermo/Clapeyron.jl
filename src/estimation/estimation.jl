include("estimationdata.jl")

struct ToEstimate
    params::Vector{Symbol}
    indices::Vector{Union{Integer,Tuple{Integer,Integer},Nothing}}  # if nothing, use all
    factor::Vector{Union{Float64,Nothing}}
    lower::Vector{Union{Vector{Union{Float64,Nothing}},Nothing}}
    upper::Vector{Union{Vector{Union{Float64,Nothing}},Nothing}}
    guess::Vector{Union{Vector{Union{Float64,Nothing}},Nothing}}  # if nothing, use current
    symmetric::Vector{Bool}
    cross_assoc::Vector{Bool}
    recombine::Vector{Bool}
end

"""
    ToEstimate
    ToEstimate(params_dict)
## Input parameters: A dictionary with the following potential entries
- `params`: The name of the parameter being fitted (`Symbol`)
- `indices`: The index of the parameter being fitted (`Integer` or `Tuple{Integer,Integer}`)
- `factor`: Factor to multiply parameter being fitted to have it in the correct units (`Float64`)
- `symmetric`: For `PairParam`, if the parameter is symmetric or asymmetric (`Bool`)
- `cross_assoc`: For `AssocParam`, if the parameter is for cross-association (`Bool`)
- `recombine`: For `PairParam`, if the combining rules must be applied for unlike interactions (`Bool`)
- `lower`: Lower bound for the parameter (`Float64`)
- `upper`: Upper bound for the parameter (`Float64`)
- `guess`: Initial guess for the parameter (`Float64`)
## Output: 
A `ToEstimate` struct
## Description
Turns the input parameter dictionary into a `ToEstimate` struct to be used within the parameter estimation.
"""
function ToEstimate(params_dict::Vector{Dict{Symbol,Any}})
    params = Vector{Symbol}(undef,0)
    indices = Vector{Union{Integer,Tuple{Integer,Integer},Nothing}}(nothing,0)
    factor = Vector{Union{Float64,Nothing}}(nothing,0)
    lower = Vector{Union{Vector{Union{Float64,Nothing}},Nothing}}(nothing,0)
    upper = Vector{Union{Vector{Union{Float64,Nothing}},Nothing}}(nothing,0)
    guess = Vector{Union{Vector{Union{Float64,Nothing}},Nothing}}(nothing,0)
    sym = Vector{Bool}(undef,0)
    recombine = Vector{Bool}(undef,0)
    cross_assoc = Vector{Bool}(undef,0)
    for dict in params_dict
        push!(params, dict[:param])
        push!(indices, get(dict, :indices, (1,1)))
        push!(factor, get(dict, :factor, 1.))
        lower_ = get(dict, :lower, nothing)
        push!(lower, typeof(lower_) <: AbstractFloat ? [lower_] : lower_)
        upper_ = get(dict, :upper, nothing)
        push!(upper, typeof(upper_) <: AbstractFloat ? [upper_] : upper_)
        guess_ = get(dict, :guess, nothing)
        push!(guess, typeof(guess_) <: AbstractFloat ? [guess_] : guess_)
        _sym = get(dict,:symmetric,true)
        push!(sym,_sym)
        _recombine = get(dict,:recombine,false)
        push!(recombine,_recombine)
        _cross_assoc = get(dict,:cross_assoc,false)
        push!(cross_assoc,_cross_assoc)
    end
    return ToEstimate(params, indices, factor, lower, upper, guess, sym, cross_assoc, recombine)
end

export Estimation
# Mutable for now to make it easy to just replace the model
"""
    Estimation
    Estimation(model,toestimate,filepaths,ignorefield)
## Input parameters:
- ` model`: The initial model containing the species we wish to parameterise
- `toestimate`: The dictionary of parameters being fitted
- `filepaths` or `filepaths_weights`: The location of the data files used to fit. Can also contain the weights of each dataset
- `ignorefield`: Specify which EoSModel fields to ignore in the main model
## Output: 
Estimator object which contains the following:
- `model`: The model whose parameters will be varied
- `initial_model`: The initial model before parameterisation
- `toestimate`: ToEstimate struct which contains all the information on the parameters
- `data`: Vector of `EstimationData` structs where all the information on the data is stored
- `ignorefield`: Vector of fields to ignore in the parameter estimation
The following objects are also output:
- `objective`: The objective function which is used to fit the parameters
- `x0`: Initial guesses for the parameters
- `upper`: Upper bounds for the parameters
- `lower`: Lower bounds for the parameters
## Description
Produces the estimator and other useful objects used within parameter estimation
"""
mutable struct Estimation{T<:EoSModel}
    model::T
    initial_model::T
    toestimate::ToEstimate
    data::Vector{EstimationData}
    ignorefield::Union{Nothing,Vector{Symbol}}
end

function Base.show(io::IO, mime::MIME"text/plain", estimation::Estimation)
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

function Base.show(io::IO, estimation::Estimation)
    print(io, typeof(estimation))
end

function Estimation(model::EoSModel, toestimate::Vector{Dict{Symbol,Any}}, filepaths::Union{Array{String},Array{Tuple{Float64, String}}}, ignorefield::Vector{Symbol})
    estimation = Estimation(model, deepcopy(model), ToEstimate(toestimate), EstimationData(filepaths),ignorefield)
    
    nparams = length(estimation.toestimate.params)

    objective(x) = objective_function(estimation,x)     
    
    x0 = [estimation.toestimate.guess[i][1] for i ∈ 1:nparams]
    upper = [estimation.toestimate.upper[i][1] for i ∈ 1:nparams]
    lower = [estimation.toestimate.lower[i][1] for i ∈ 1:nparams]
    return estimation, objective, x0, upper, lower
end

function Estimation(model::EoSModel, toestimate::Vector{Dict{Symbol,Any}}, filepaths::Union{Array{String},Array{Tuple{Float64, String}}})
    estimation = Estimation(model, deepcopy(model), ToEstimate(toestimate), EstimationData(filepaths),Symbol[])
    
    nparams = length(estimation.toestimate.params)

    objective(x) = objective_function(estimation,x)
    
    x0 = [estimation.toestimate.guess[i][1] for i ∈ 1:nparams]
    upper = [estimation.toestimate.upper[i][1] for i ∈ 1:nparams]
    lower = [estimation.toestimate.lower[i][1] for i ∈ 1:nparams]
    return estimation, objective, x0, upper, lower
end

function reload_data(estimation::Estimation)
    estimationdata = EstimationData(estimation.filepaths)
    empty!(estimation.data)
    for i in 1:length(estimation.filepaths)
        push!(estimation.data, estimationdata[i])
    end
end

export update_estimation!
function update_estimation!(
        estimation::Estimation,
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

function update_estimation!(estimation::Estimation, model::EoSModel)
    estimation.model = model
end

export return_model
"""
    return_model
    return_model(estimation,model,params)
## Input parameters:
- `estimation`: The estimator object
- `model`: The model whose parameters we are varying
- `params`: The new parameters which we want to change
## Output: 
- `model`: The new model with the updated parameters
## Description
Based on the parameters provided and the estimator, a new model is produced from the input.
"""
function return_model(
        estimation::Estimation,
        model::EoSModel,
        values) 
    params = estimation.toestimate.params
    factor = estimation.toestimate.factor
    sym = estimation.toestimate.symmetric
    cross_assoc = estimation.toestimate.cross_assoc
    idx = estimation.toestimate.indices
    recombine = estimation.toestimate.recombine
    model = deepcopy(model)
    if isdefined(model,:params)
        for (i, param) in enumerate(params)
            f = factor[i]
            id = idx[i]
            recomb = recombine[i]
            if isdefined(model.params,param)
                current_param = getfield(model.params, param)
                if typeof(current_param) <: SingleParameter
                    current_param[id[1]] = values[i]*f
                end
                if typeof(current_param) <: PairParam
                    current_param[id[1],id[2],sym[i]] = values[i]*f
                    if (id[1]==id[2]) & recomb
                        current_param.ismissingvalues[id[1],:] .= true
                        current_param.ismissingvalues[:,id[1]] .= true
                    elseif id[1]!=id[2]
                        current_param.ismissingvalues[id[1],id[2]] = false
                        current_param.ismissingvalues[id[2],id[1]] = false
                    end
                end
                if typeof(current_param) <: AssocParam
                    current_param.values.values[id[1]] = values[i]*f
                    if cross_assoc[i]
                        current_param.values.values[id[1]+1] = values[i]*f                
                    end
                end
            end
        end
    end
    for i ∈ fieldnames(typeof(model))
        if (typeof(getfield(model,i)) <: EoSModel) & !(i in estimation.ignorefield)
            return_model!(estimation,getfield(model,i),values)
        end
    end
    recombine!(model)
    return model
end

function return_model!(
    estimation::Estimation,
    model::EoSModel,
    values) 
    params = estimation.toestimate.params
    factor = estimation.toestimate.factor
    sym = estimation.toestimate.symmetric
    cross_assoc = estimation.toestimate.cross_assoc
    idx = estimation.toestimate.indices
    recombine = estimation.toestimate.recombine
    if isdefined(model,:params)
        for (i, param) in enumerate(params)
            f = factor[i]
            id = idx[i]
            recomb = recombine[i]
            if isdefined(model.params,param)
                current_param = getfield(model.params, param)
                if typeof(current_param) <: SingleParameter
                    current_param[id[1]] = values[i]*f
                end
                if typeof(current_param) <: PairParam
                    current_param[id[1],id[2],sym[i]] = values[i]*f
                    if (id[1]==id[2]) & recomb
                        current_param.ismissingvalues[id[1],:] .= true
                        current_param.ismissingvalues[:,id[1]] .= true
                    elseif id[1]!=id[2]
                        current_param.ismissingvalues[id[1],id[2]] = false
                        current_param.ismissingvalues[id[2],id[1]] = false
                    end
                end
                if typeof(current_param) <: AssocParam
                    current_param.values.values[id[1]] = values[i]*f
                    if cross_assoc[i]
                        current_param.values.values[id[1]+1] = values[i]*f                
                    end
                end
            end
        end
    end
    for i ∈ fieldnames(typeof(model))
        if (typeof(getfield(model,i)) <: EoSModel) & !(i in estimation.ignorefield)
            return_model!(estimation,getfield(model,i),values)
        end
    end
    recombine!(model)
end

"""
    objective_function
    objective_function(estimation,params)
## Input parameters:
- `estimation`: The estimator object
- `params`: The new parameters which we want to evaluate the objective function for
## Output: The relate root mean square error given the data and parameters provided
## Description
The objective function used within parameter estimation.
"""
function objective_function(estimation::Estimation,guesses)
    F = 0
    model = return_model(estimation, estimation.model, guesses)
    for i ∈ 1:length(estimation.data)
        property = estimation.data[i].method
        inputs = estimation.data[i].inputs
        outputs = estimation.data[i].outputs
        weights = estimation.data[i].weights
        if isempty(inputs)
            prediction =  property(model)
        elseif length(inputs)==1
            prediction =  property.(model,inputs[1])
        elseif length(inputs)==2
            prediction =  property.(model,inputs[1],inputs[2])
        elseif length(inputs)==3
            prediction =  property.(model,inputs[1],inputs[2],inputs[3])
        else
            prediction = property.(model,inputs...)
        end

        if length(outputs)==1
            F += sum(((prediction.-outputs[1])./outputs[1]).^2)/length(outputs[1])*weights[1]
        else
            F += sum([sum([((prediction[k][j].-outputs[j][i])./outputs[j][k]).^2 for j in 1:length(prediction[k])])*weights[1] for k in 1:length(prediction)])/length(outputs[1])
        end
    end
    if isnan(F)
        return 1e4
    else
        return F
    end
    return_model!(estimation,estimation.model,guesses)
end