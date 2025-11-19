function init_preferred_method(method::typeof(bubble_pressure),model::RestrictedEquilibriaModel,kwargs)
    return ActivityBubblePressure(;kwargs...)
end