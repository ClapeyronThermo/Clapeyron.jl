function excess_gibbs_free_energy(model::ActivityModel,p,T,z)
    lnγ = log.(activity_coefficient(model,p,T,z))
    return sum(z[i]*R̄*T*lnγ[i] for i ∈ @comps)
end

function bubble_pressure(model::ActivityModel,T,x)
    sat = sat_pure.(model.puremodel,T)
    p_sat = [tup[1] for tup in sat]
    γ     = activity_coefficient(model,1e-4,T,x)
    p     = sum(x.*γ.*p_sat)
    y     = x.*γ.*p_sat ./ p
    return (p,y)
end
