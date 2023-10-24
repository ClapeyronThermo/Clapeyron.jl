function saturation_p(model::EoSModel,T)
    sat = saturation_pressure(model,T)
    return sat[1]
end

function saturation_rhol(model::EoSModel,T)
    sat = saturation_pressure(model,T)
    return 1/sat[2]
end

function crit_temp(model::EoSModel)
    crit = crit_pure(model)
    return crit[1]
end

function bubble_temp(model::EoSModel,p,x)
    x = Fractions.FractionVector(x)
    bub = bubble_temperature(model,p,x)
    return bub[1]
end

function bubble_comp(model::EoSModel,p,x)
    x = Fractions.FractionVector(x)
    bub = bubble_temperature(model,p,x)
    return bub[4][1]
end