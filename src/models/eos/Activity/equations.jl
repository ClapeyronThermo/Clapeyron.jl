function excess_gibbs_free_energy(model::ActivityModel,p,T,z)
    γ = activity_coefficient(model,p,T,z)
    return sum(z[i]*R̄*T*log(γ[i]) for i ∈ @comps)
end

function eos(model::ActivityModel,V,T,z)
    Σz = sum(z)
    lnΣz = log(Σz)
    g_E = excess_gibbs_free_energy(model,V,T,z)
    pures = model.puremodel
    p = sum(z[i]*pressure(pures[i],V,T) for i ∈ @comps)/Σz
    g_ideal = sum(z[i]*R̄*T*(log(z[i])-lnΣz) for i ∈ @comps)
    g_pure = sum(z[i]*VT_gibbs_free_energy(pures[i],V,T) for i ∈ @comps)
    return g_E+g_ideal+g_pure-p*V
end

function eos_res(model::ActivityModel,V,T,z)
    Σz = sum(z)
    g_E = excess_gibbs_free_energy(model,V,T,z)
    pures = model.puremodel
    g_pure_res = sum(VT_chemical_potential_res(pures[i],V,T)[1]*z[i] for i ∈ @comps)
    p = sum(z[i]*pressure(pures[i],V,T) for i ∈ @comps)/Σz
    p_res = p - Σz*R̄*T/V
    return g_E+g_pure_res-p_res*V
end

function bubble_pressure(model::ActivityModel,T,x)
    sat = saturation_pressure.(model.puremodel,T)
    p_sat = [tup[1] for tup in sat]
    γ     = activity_coefficient(model,1e-4,T,x)
    p     = sum(x.*γ.*p_sat)
    y     = x.*γ.*p_sat ./ p
    return (p,y)
end

function bubble_temperature(model::ActivityModel,p,x;T0=nothing)
    f(z) = Obj_bubble_temperature(model,z,p,x)
    if T0===nothing
        pure = model.puremodel
        sat = saturation_temperature.(pure,p)
        Ti   = zero(x)
        for i ∈ 1:length(x)
            if isnan(sat[i][1])
                Tc,pc,vc = crit_pure(pure[i])
                g(x) = p-pressure(pure[i],vc,x,[1.])
                Ti[i] = Roots.find_zero(g,(Tc))
            else
                Ti[i] = sat[i][1]
            end
        end
        T = Roots.find_zero(f,(minimum(Ti)*0.9,maximum(Ti)*1.1))
    else
        T = Roots.find_zero(f,T0)
    end
    p,y = bubble_pressure(model,T,x)
    return (T,y)
end

function Obj_bubble_temperature(model::ActivityModel,T,p,x)
    sat = saturation_pressure.(model.puremodel,T)
    p_sat = [tup[1] for tup in sat]
    γ     = activity_coefficient(model,1e-4,T,x)
    y     = x.*γ.*p_sat ./ p
    return sum(y)-1
end