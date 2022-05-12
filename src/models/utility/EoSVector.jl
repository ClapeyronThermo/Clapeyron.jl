struct EoSVector{T} <: EoSModel
    components::Vector{String}
    mix::T
    pure::Vector{T}
end

function EoSVector(model::EoSModel)
    pure = split_model(model)
    components = model.components
    return EoSVector(components,model,pure)
end

Base.getindex(x::EoSVector,I) = x.pure[I]
Base.length(x::EoSVector) = length(x.pure)
Base.eltype(x::EoSVector{T}) where T = T 
Base.broadcastable(x::EoSVector) = x.pure

function init_puremodel(model::Type{<:EoSModel},components,userlocations,verbose)
    verbose && @info("""Now creating ideal model:
    $idealmodel""")
    return EoSVector(model(components;userlocations,verbose))
end

function init_puremodel(model::EoSModel,components,userlocations,verbose)
   return EoSVector(model)
end

eos(model::EoSVector,V,T,z=SA[1.0]) = eos(model.mix,V,T,z)
eos_res(model::EoSVector,V,T,z=SA[1.0]) = eos_res(model.mix,V,T,z)
idealmodel(model::EoSVector) = idealmodel(model.mix)

