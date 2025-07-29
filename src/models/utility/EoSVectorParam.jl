struct EoSVectorParam{T} <: EoSModel
    components::Vector{String}
    model::T
    pure::Vector{T}
    reference_state::Union{Nothing,ReferenceState}
end

function EoSVectorParam(model::EoSModel,components = model.components)
    pure = split_pure_model(model,1:length(components))
    if has_reference_state(model)
        ref = nothing
    else
        ref = ReferenceState()
    end
    return EoSVectorParam(components,model,pure,ref)
end

saturation_model(model::EoSVectorParam) = saturation_model(model.model)
a_res(model::EoSVectorParam,V,T,z) = a_res(model.model,V,T,z)
a_res(model::EoSVectorParam,V,T,z,data) = a_res(model.model,V,T,z,data)
Rgas(model::EoSVectorParam) = Rgas(model.model)
molecular_weight(model::EoSVectorParam,z) = molecular_weight(model.model,z)
function eos_impl(model::EoSVectorParam,V,T,z)
    return Rgas(model)*sum(z)*T*a_eos(model.model,V,T,z) + reference_state_eval(model,V,T,z)
end

Base.getindex(x::EoSVectorParam,I) = x.pure[I]
Base.length(x::EoSVectorParam) = length(x.pure)
Base.eltype(x::EoSVectorParam{T}) where T = T
Base.broadcastable(x::EoSVectorParam) = x.pure

function init_puremodel(model,components,userlocations,verbose)
    verbose && @info("Now creating pure model: " * info_color(string(model)))
    _model = init_model(model,components,userlocations,verbose)
    return init_puremodel(_model,components,userlocations,verbose)
end

function init_puremodel(model::EoSModel,components,userlocations,verbose)
    _components = format_components(components)
    pure = split_model(model,1:length(_components))
    if has_reference_state(model)
        ref = nothing
    else
        ref = ReferenceState()
    end
    return EoSVectorParam(_components,model,pure,ref)
end

function recombine_impl!(model::EoSVectorParam)
    recombine!(model.model)
    if is_splittable(model.model)
        model.pure .= split_pure_model(model.model)
    else
        model.pure .= model.model
    end
    return model
end

function PT_property(model::EoSVectorParam,p,T,z,phase,threaded,vol0,f::F,::Val{UseP}) where {F,UseP}
    v = volume(model.model,p,T,z;phase,threaded,vol0)
    if UseP
        return f(model,v,T,z,p)
    else
        return f(model,v,T,z)
    end
end

function reference_state(model::EoSVectorParam)
    original_ref = reference_state(model.model)
    if original_ref == nothing
        return model.reference_state
    else
        return original_ref
    end
end

function volume_impl(model::EoSVectorParam,p,T,z,phase,threaded,vol0)
    return volume_impl(model.model,p,T,z,phase,threaded,vol0)
end

function promote_model(::Type{T},model::EoSVectorParam) where T <: Number
    return EoSVectorParam(model.components,promote_model(T,model.model),promote_model(T,model.pure),model.reference_state)
end