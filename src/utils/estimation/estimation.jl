include("estimationdata.jl")

export Estimation
# Mutable for now to make it easy to just replace the model
mutable struct Estimation{T<:EoSModel}
    model::T
    initial_model::T
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
end

function Base.show(io::IO, estimation::Estimation)
    print(io, typeof(estimation))
end

function Estimation(model::EoSModel, filepaths::Array{String})
    return Estimation(model, deepcopy(model), filepaths, EstimationData(filepaths))
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
        params::Vector{Symbol},
        values::Vector{T} where {T<:Any}) 
    model = deepcopy(estimation.model)
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
    return model
end

