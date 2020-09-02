function a_res(model::PcSaftFamily, conditions)
    return a_hc(model, conditions) + a_disp(model, conditions)
end

function a_hc(model::PcSaftFamily, conditions)
    x = conditions.compositions
    m = model.parameters.segments
    m̄ = sum(x[i]*m[i] for i in model.components)
    return m̄*a_hs(model, conditions) - sum(x[i]*(m[i]-1)*log(g_hsij(model, conditions, i, i)) for i in model.components)
end

function a_disp(model::PcSaftFamily, conditions)
    v = conditions.volume
    pai = 3.1415
    return -2*pai/v*I(model, conditions, 1)*m2ϵσ3(model, conditions, 1) - pai/v*C1(model, conditions)*I(model, conditions, 2)*m2ϵσ3(model, conditions, 2)
end

function d(model::PcSaftFamily, conditions, component)
    T = conditions.temperature
    ϵ = model.parameters.epsilons[component]
    σ = model.parameters.sigmas[component]
    return σ * (1 - 0.12exp(-3ϵ/T))
end
    
function ζn(model::PcSaftFamily, conditions, n)
    v = conditions.volume
    x = conditions.compositions
    m = model.parameters.segments
    pai = 3.1415
    return pai/6/v * sum(x[i]*m[i]*d(model, conditions, i)^n for i in model.components)
end

function g_hsij(model::PcSaftFamily, conditions, i, j)    
    di = d(model, conditions, i)
    dj = d(model, conditions, j)
    ζ2 = ζn(model, conditions, 2)
    ζ3 = ζn(model, conditions, 3)
    return 1/(1-ζ3) + di*dj/(di+dj)*3ζ2/(1-ζ3)^2 + (di*dj/(di+dj))^2*2ζ2^2/(1-ζ3)^3
end

function a_hs(model::PcSaftFamily, conditions)
    ζ0 = ζn(model, conditions, 0)
    ζ1 = ζn(model, conditions, 1)
    ζ2 = ζn(model, conditions, 2)
    ζ3 = ζn(model, conditions, 3)
    return 1/ζ0 * (3ζ1*ζ2/(1-ζ3) + ζ2^3/(ζ3*(1-ζ3)^2) + (ζ2^3/ζ3^2-ζ0)*log(1-ζ3))
end

function C1(model::PcSaftFamily, conditions)
    x = conditions.compositions
    η = ζn(model, conditions, 3)
    m = model.parameters.segments
    m̄ = sum(x[i]*m[i] for i in model.components)
    return (1 + m̄*(8η-2η^2)/(1-η)^4 + (1-m̄)*(20η-27η^2+12η^3-2η^4)/((1-η)*(2-η))^2)^-1
end

function m2ϵσ3(model::PcSaftFamily, conditions, ϵ_power = 1)
    T = conditions.temperature
    x = conditions.compositions
    m = model.parameters.segments
    σ = model.parameters.sigmas
    ϵ = model.parameters.epsilons
    k = model.parameters.ks
    return sum(x[i]*x[j]*m[i]*m[j] * (sqrt(ϵ[i]*ϵ[j])*(1-k[(i,j)])/T)^ϵ_power * (0.5*(σ[i]+σ[j]))^3 for i in model.components, j in model.components)
end

function I(model::PcSaftFamily, conditions, n)
    x = conditions.compositions
    m = model.parameters.segments
    m̄ = sum(x[i]*m[i] for i in model.components)
    η = ζn(model, conditions, 3)
    if n == 1
        corr = [0.9105631445 -0.3084016918 -0.0906148351;
                0.6361281449 0.1860531159 0.4527842806;
                2.6861347891 -2.5030047259 0.5962700728;
                -26.547362491 21.419793629 -1.7241829131;
                97.759208784 -65.255885330 -4.1302112531;
                -159.59154087 83.318680481 13.776631870;
                91.297774084 -33.746922930 -8.6728470368]
    elseif n == 2
        corr = [0.7240946941 -0.5755498075 0.0976883116;
                2.2382791861 0.6995095521 -0.2557574982;
                -4.0025849485 3.8925673390 -9.1558561530;
                -21.003576815 -17.215471648 20.642075974;
                26.855641363 192.67226447 -38.804430052;
                206.55133841 -161.82646165 93.626774077;
                -355.60235612 -165.20769346 -29.666905585]
    end
    return sum((corr[i,1] + (m̄-1)/m̄*corr[i,2] + (m̄-1)/m̄*(m̄-2)/m̄*corr[i,3]) * η^i for i = 1:6)
end
