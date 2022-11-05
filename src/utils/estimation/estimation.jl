include("estimationdata.jl")

struct ToEstimate
    params::Vector{Symbol}
    indices::Vector{Union{Vector{Integer},Nothing}}  # if nothing, use all
    factor::Vector{Union{Float64,Nothing}}
    lower::Vector{Union{Vector{Union{Float64,Nothing}},Nothing}}
    upper::Vector{Union{Vector{Union{Float64,Nothing}},Nothing}}
    guess::Vector{Union{Vector{Union{Float64,Nothing}},Nothing}}  # if nothing, use current
    symmetric::Vector{Bool}
end

function ToEstimate(params_dict::Vector{Dict{Symbol,Any}})
    params = Vector{Symbol}(undef,0)
    indices = Vector{Union{Vector{Integer},Nothing}}(nothing,0)
    factor = Vector{Union{Float64,Nothing}}(nothing,0)
    lower = Vector{Union{Vector{Union{Float64,Nothing}},Nothing}}(nothing,0)
    upper = Vector{Union{Vector{Union{Float64,Nothing}},Nothing}}(nothing,0)
    guess = Vector{Union{Vector{Union{Float64,Nothing}},Nothing}}(nothing,0)
    sym = Vector{Bool}(undef,0)
    for dict in params_dict
        push!(params, dict[:param])
        push!(indices, get(dict, :indices, nothing))
        push!(factor, get(dict, :factor, 1.))
        lower_ = get(dict, :lower, nothing)
        push!(lower, typeof(lower_) <: AbstractFloat ? [lower_] : lower_)
        upper_ = get(dict, :upper, nothing)
        push!(upper, typeof(upper_) <: AbstractFloat ? [upper_] : upper_)
        guess_ = get(dict, :guess, nothing)
        push!(guess, typeof(guess_) <: AbstractFloat ? [guess_] : guess_)
        _sym = get(dict,:symmetric,true)
        push!(sym,_sym)
    end
    return ToEstimate(params, indices, factor, lower, upper, guess, sym)
end

export Estimation
# Mutable for now to make it easy to just replace the model
mutable struct Estimation{T<:EoSModel}
    model::T
    initial_model::T
    toestimate::ToEstimate
    filepaths::Array{String}
    data::Vector{EstimationData}
end

function Base.show(io::IO, mime::MIME"text/plain", estimation::Estimation)
    print(io, typeof(estimation))
    println(io, " with data for:")
    firstloop = true
    for data in estimation.data
        !firstloop && println(io, "")
        print(io, "  :" * String(data.method))
        firstloop = false
    end
    println(io, "\n to estimate:")
    firstloop = true
    for (param, indices) in zip(
            estimation.toestimate.params, estimation.toestimate.indices)
        !firstloop && println(io, "")
        print(io, "  :" * String(param))
        if !(indices === nothing)
            print(io, " with indices => " * "[" * 
                  join(indices,",") * "]")
        end
        firstloop = false
    end
end

function Base.show(io::IO, estimation::Estimation)
    print(io, typeof(estimation))
end

function Estimation(model::EoSModel, toestimate::Vector{Dict{Symbol,Any}}, filepaths::Array{String})
    return Estimation(model, deepcopy(model), ToEstimate(toestimate), filepaths, EstimationData(filepaths))
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
function return_model(
        estimation::Estimation,
        values::Vector{T} where {T<:Any}) 
    params = estimation.toestimate.params
    factor = estimation.toestimate.factor
    sym = estimation.toestimate.symmetric
    model = deepcopy(estimation.model)
    for (i, param) in enumerate(params)
        f = factor[i]
        current_param = getfield(model.params, param)
        if typeof(current_param) <: SingleParameter
            for (j, value) in enumerate(values[i])
                current_param[j] = value*f
            end
        end
        if typeof(current_param) <: PairParam
            for j = 1:length(model.components)
                for k = 1:length(model.components)
                    current_param[j,k,sym[i]] = values[i][j,k]*f
                end
            end
        end
        if typeof(current_param) <: AssocParam
            if typeof(values[i]) <: Compressed4DMatrix
                for (j, value) in enumerate(values[i].values)
                    current_param.values.values[j] = value*f
                end
            else
                for (j, value) in enumerate(values[i])
                    current_param.values.values[j] = value*f
                end
            end
        end
    end
    return model
end

export optimize!

function optimize!(estimation::Estimation,Method=:random_search,MaxSteps=20000,PopulationSize=500)
    f(x) = obj_fun(estimation,x)
    nparams = length(estimation.toestimate.params)
    bounds = [(estimation.toestimate.lower[i][1],estimation.toestimate.upper[i][1]) for i ∈ 1:nparams]
    x0 = [estimation.toestimate.guess[i][1] for i ∈ 1:nparams]
    x=BlackBoxOptim.bboptimize(f; SearchRange = bounds, NumDimensions=nparams, Method=Method, MaxSteps=MaxSteps, PopulationSize=PopulationSize);
    model = return_model(estimation, BlackBoxOptim.best_candidate(x))
    update_estimation!(estimation,model)
end

function obj_fun(estimation::Estimation,guesses)
    F = 0
    model = return_model(estimation, guesses)
    for i ∈ 1:length(estimation.data)
        property = getfield(Clapeyron,estimation.data[i].method)
        inputs = estimation.data[i].inputs
        output = estimation.data[i].outputs[1]
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
        F += sum(abs.((prediction.-output)./output))/length(output)
    end
    return F
end

toestimate = [
    Dict(
        :param => :epsilon,
        :indices => [1],
        :lower => [3.7],
        :upper => 5.0,
        :guess => 3.0
    ),
    Dict(
        :param => :sigma,
        :lower => 3.3,
        :upper => 3.8,
        :guess => 3.5
    ),
    Dict(
        :param => :lambda_r,
        :lower => 12.0,
        :upper => 16.0,
        :guess => 16.0
    )
] 
