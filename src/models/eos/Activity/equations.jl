function excess_gibbs_free_energy(model::ActivityModel,p,T,z)
    lnγ = log.(activity_coefficient(model,p,T,z))
    return sum(z[i]*R̄*T*lnγ[i] for i ∈ @comps)
end

function eos(model::ActivityModel,V,T,z)
    x = z./sum(z)
    g_E = excess_gibbs_free_energy(model,V,T,z)
    g_pure = VT_gibbs_free_energy.(model.puremodel,V,T)
    p      = pressure.(model.puremodel,V,T)
    g_ideal = sum(z[i]*R̄*T*log(x[i]) for i ∈ @comps)
    return g_E+g_ideal+sum(z[i]*g_pure[i] for i ∈ @comps)-sum(x[i]*p[i] for i ∈ @comps)*V
end

function eos_res(model::ActivityModel,V,T,z)
    x = z./sum(z)
    g_E = excess_gibbs_free_energy(model,V,T,z)
    g_pure_res = VT_chemical_potential_res.(model.puremodel,V,T)
    p_res      = pressure.(model.puremodel,V,T).-sum(z)*R̄*T/V
    return g_E+sum(z[i]*g_pure_res[i][1] for i ∈ @comps)-sum(z[i]*p_res[i] for i ∈ @comps)*V
end

function bubble_pressure(model::ActivityModel,T,x)
    sat = sat_pure.(model.puremodel,T)
    p_sat = [tup[1] for tup in sat]
    γ     = activity_coefficient(model,1e-4,T,x)
    p     = sum(x.*γ.*p_sat)
    y     = x.*γ.*p_sat ./ p
    return (p,y)
end
