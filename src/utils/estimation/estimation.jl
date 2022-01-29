export Estimation
# Mutable for now to make it easy to just replace the model
mutable struct Estimation{T<:EoSModel, U<:EoSParam}
    model::T
    initial_param::U
    # data::Vector{EstimationData}
end

function Estimation(model::EoSModel)
    return Estimation(model, deepcopy(model.params))
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

