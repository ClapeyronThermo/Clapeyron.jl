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
    if is_splittable(_model)
        pure = split_model(_model)
    else
        pure = fill(_model,length(components))
    end
        return EoSVectorParam(components,_model,pure)
end



function init_puremodel(model::EoSModel,components,userlocations,verbose)
   return EoSVectorParam(model)
end


