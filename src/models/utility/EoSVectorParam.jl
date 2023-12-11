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
    verbose && @info("Now creating pure model: " * info_color(string(model)))
    _model = init_model(model,components,userlocations,verbose)
    return init_puremodel(_model,components,userlocations,verbose)
end

function init_puremodel(model::EoSModel,components,userlocations,verbose)
    _components = format_components(components)
    pure = split_model(model,1:length(_components))
    return EoSVectorParam(_components,model,pure)
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

function saturation_pressure(model::EoSVectorParam,T::Real,method::SaturationMethod)
    return saturation_pressure(model.model,T,method)
end

function saturation_temperature(model::EoSVectorParam,T::Real,method::SaturationMethod)
    return saturation_temperature(model.model,T,method)
end

function init_preferred_method(method::typeof(saturation_pressure),model::EoSVectorParam,kwargs)
    return init_preferred_method(method,model.model,kwargs)
end

function init_preferred_method(method::typeof(saturation_temperature),model::EoSVectorParam,kwargs)
    return init_preferred_method(method,model.model,kwargs)
end
