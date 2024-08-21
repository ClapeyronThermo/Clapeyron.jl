function equilibrate(model::ReactiveAqModel,p,T,n0;z0=nothing)
    ν = stoichiometric_coefficient(model) # Exists
    Keq = ideal_Keq(model,T,n0,ν)
    println(Keq)
    if isnothing(z0)
        z0 = x0_equilibrium_conditions(model,p,T,n0)
    end
    f!(F,z) = equilibrium_conditions(model,F,p,T,n0,exp10.(z),ν,Keq)
    sol = Solvers.nlsolve(f!,log10.(z0),TrustRegion(Newton(), NWI()),NEqOptions(),ForwardDiff.Chunk{1}()) # Exists
    ξ = exp10.(Solvers.x_sol(sol))
    n = n0+ν*ξ
    return n
end

function pH(model::ReactiveAqModel,p,T,n0;z0=nothing)
    n = equilibrate(model, p, T, n0; z0 = z0)
    a = aqueous_activity(model.eosmodel,p,T,n)
    hydronium_id = find_hydronium_index(model)
    hydroxide_id = find_hydroxide_index(model)

    idx_w = find_water_indx(model)

    zref = ones(length(model.components)).*1e-30
    zref[idx_w] = 1.
    zref[hydronium_id] = 0.01801528*1.
    zref[hydroxide_id] = 0.01801528*1.

    zref ./= sum(zref)

    μref = chemical_potential(model.eosmodel,p,T,zref)
    μmix = chemical_potential(model.eosmodel,p,T,n)

    a = exp.((μmix .- μref)./R̄./T)

    return -log10(a[hydronium_id])
end

function equilibrium_conditions(model::ReactiveAqModel,F,p,T,n0,ξ,ν,Keq)
    n = n0+ν*ξ

    nreact = length(model.reactions)

    μmix = chemical_potential(model.eosmodel,p,T,n)

    for i in 1:nreact
        idx_w = find_water_indx(model)

        zref = ones(length(model.components)).*1e-30
        zref[idx_w] = 1.
        zref[ν[:,i] .!= 0 .&& model.eosmodel.charge .!= 0] .= 0.01801528*1.

        zref ./= sum(zref)

        μref = chemical_potential(model.eosmodel,p,T,zref)

        a = exp.((μmix .- μref)./R̄./T)
        F[i] = sum(log.(a).*ν[:,i]) - log(Keq[i])
    end
    return F
end