function pH(model::ReactiveAqModel,p,T,n0;z0=nothing)
    n = equilibrate(model, p, T, n0; z0 = z0)

    a = aqueous_activity(model.eosmodel,p,T,n)

    hydronium_id = find_hydronium_index(model)

    return -log10(a[hydronium_id])
end

function equilibrate(model::ReactiveAqModel,p,T,n0;z0=nothing)
    T0 = 298.15
    ν = stoichiometric_coefficient(model)
    ΔGf0 = model.params.ΔHf.values
    ΔHf = model.params.ΔHf.values
    ΔGf = ΔGf0/T0+ΔHf*(1/T-1/T0)
    ΔrG = sum(ΔGf.*ν,dims=1)
    Keq = exp.(-ΔrG/Rgas(model.eosmodel))

    # println(Keq)

    if isnothing(z0)
        z0 = x0_equilibrium_conditions(model,p,T,n0)
    end

    f!(F,z) = equilibrium_conditions(model,F,p,T,n0,exp10.(z),ν,Keq)
    sol = Solvers.nlsolve(f!,log10.(z0),TrustRegion(Newton(), NWI()),NEqOptions(),ForwardDiff.Chunk{1}()) # Exists
    ξ = exp10.(Solvers.x_sol(sol))
    n = n0+ν*ξ
    return n
end

function equilibrium_conditions(model::ReactiveAqModel,F,p,T,n0,ξ,ν,Keq)
    n = n0+ν*ξ
    # println(n)

    a = aqueous_activity(model.eosmodel,p,T,n)

    F[1:end] = sum(log.(a).*ν,dims=1) - log.(Keq)
    # println(F)
    return F
end

function x0_equilibrium_conditions(model::ReactiveAqModel,p,T,n0)
    
    z0 = zeros(length(model.components))
    z0[1] = log10(p)
    return z0
end

export equilibrate