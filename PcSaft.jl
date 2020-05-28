#= function a_ideal(model::PcSaftFamily, conditions) =#
#=     0 =#
#= end =#
#= function a_mono(model::PcSaftFamily, conditions) =#
#=     0 =#
#= end =#
#= function a_chain(model::PcSaftFamily, conditions) =#
#=     0 =#
#= end =#
#= function a_assoc(model::PcSaftFamily, conditions) =#
#=     0 =#
#= end =#

function d(model::PcSaftFamily, conditions, component)
    T = conditions.temperature
    epsilon = model.parameters.epsilons[component]
    sigma = model.parameters.sigmas[component]
    return sigma * (1 - 0.12exp(-3epsilon/T))
end
    
function zetan(model::PcSaftFamily, conditions, n)
    v = conditions.volume
    components = conditions.components
    segments = model.parameters.segments
    return pi/6/v*sum(components[i]*segments[i]*d(model, conditions, i)^n for i in keys(components))
end

function ghsij(model::PcSaftFamily, conditons, i, j)    
    di = d(model, conditons, i)
    dj = d(model, conditons, j)
    zeta2 = zetan(model, conditons, 2)
    zeta3 = zetan(model, conditons, 3)
    return 1/(1-zeta3) + di*dj/(di+dj)*3zeta2/(1-zeta3)^2 + (di*dj/(di+dj))^2*2zeta2^2/(1-zeta3)^3
end

function ahs(model::PcSaftFamily, conditons)
    zeta0 = zetan(model, conditons, 0)
    zeta1 = zetan(model, conditons, 1)
    zeta2 = zetan(model, conditons, 2)
    zeta3 = zetan(model, conditons, 3)
    return 1/zeta0*( 3zeta1*zeta2/(1-zeta3) + zeta2^3/(zeta3*(1-zeta3)^2) + (zeta2^3/zeta3^2-zeta0)*log(1-zeta3) )
end

function ahc(model::PcSaftFamily, conditons)
    components = conditions.components
    segments = model.parameters.segments
    mbar = sum(components[i]*segments[i] for i in keys(components))
    return mbar*ahs(model,conditions) - sum(components[i]*(segments[i]-1)*log(ghsij(model, conditions, i, i)) for i in keys(components))
end

function adisp(model::PcSaftFamily, conditons)
    return 0
end
