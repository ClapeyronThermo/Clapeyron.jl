function init_preferred_method(method::typeof(dew_pressure),model::RestrictedEquilibriaModel,kwargs)
    return ActivityDewPressure(;kwargs...)
end




function init_preferred_method(method::typeof(dew_temperature),model::RestrictedEquilibriaModel,kwargs)
    return FugDewTemperature(;kwargs...)
end