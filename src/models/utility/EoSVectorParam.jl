struct EoSVectorParam{T} <: EoSModel
    components::Vector{String}
    model::T
    pure::Vector{T}
end

function EoSVectorParam(model::EoSModel)
    pure = split_model(model)
    components = model.components
    return EoSVectorParam(components,model,pure)
end

Base.getindex(x::EoSVectorParam,I) = x.pure[I]
Base.length(x::EoSVectorParam) = length(x.pure)
Base.eltype(x::EoSVectorParam{T}) where T = T
Base.broadcastable(x::EoSVectorParam) = x.pure

function init_puremodel(model,components,userlocations,verbose)
    verbose && @info("""Now creating pure model:
    $idealmodel""")
    _model = init_model(model,components,userlocations,verbose)
    return init_puremodel(_model,components,userlocations,verbose)
end

function init_puremodel(model::EoSModel,components,userlocations,verbose)
    if is_splittable(model)
        pure = split_model(model)
    else
        pure = fill(model,length(components))
    end
    return EoSVectorParam(components,model,pure)
end

function recombine_impl!(model::EoSVectorParam)
    recombine!(model.model)
    if is_splittable(model.model)
        model.pure .= split_model(model.model)
    else
        model.pure .= model.model
    end
    return model
end