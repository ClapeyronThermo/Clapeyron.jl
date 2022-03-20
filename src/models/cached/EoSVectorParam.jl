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

function init_puremodel(model::Type{<:EoSModel},components,userlocations,verbose)
    verbose && @info("""Now creating ideal model:
    $idealmodel""")
    return EoSVectorParam(model(components;userlocations,verbose))
end

function init_puremodel(model::EoSModel,components,userlocations,verbose)
   return EoSVectorParam(model)
end
