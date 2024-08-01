# function pH(model::ReactiveModel,p,T,n0;z0=nothing)
#     ν = stoichiometric_coefficient(model) # Exists
#     ΔHf = model.params.ΔHf.values
#     Sf = model.params.Sf.values
#     ΔGf = ΔHf-T*Sf
#     ΔrG = sum(ΔGf*ν,dims=1)
#     Keq = exp.(-ΔrG/Rgas(model.model)/T)

#     if isnothing(z0)
#         z0 = x0_equilibrium_conditions(model,p,T,n0)
#     end

#     f!(F,z) = equilibrium_conditions(model,F,p,T,n0,z,ν,Keq)
#     sol = Solvers.nlsolve(f!,z0,TrustRegion(Newton(), NWI()),NEqOptions(),ForwardDiff.Chunk{1}()) # Exists
#     ξ = Solvers.x_sol(sol)
#     n = n0+ν*ξ

#     a = activity(model.model,p,T,n) # Exists

#     hydronium_id = find_hydronium_index(model)

#     return -log10(a[hydronium_id])
# end

# function find_hydronium_index(model)
#     idx = findfirst(model.components.=="hydronium")
#     return idx
# end

function equilibrate(model::ReactiveModel,p,T,n0;z0=nothing)
    ν = stoichiometric_coefficient(model) # Exists
    ΔHf = model.params.ΔHf.values
    Sf = model.params.Sf.values
    ΔGf = ΔHf-T*Sf
    ΔrG = sum(ΔGf.*ν,dims=1)
    Keq = exp.(-ΔrG/Rgas(model.eosmodel)/T)

    if isnothing(z0)
        z0 = x0_equilibrium_conditions(model,p,T,n0)
    end

    f!(F,z) = equilibrium_conditions(model,F,p,T,n0,z,ν,Keq)
    sol = Solvers.nlsolve(f!,z0,TrustRegion(Newton(), NWI()),NEqOptions(),ForwardDiff.Chunk{1}()) # Exists
    ξ = Solvers.x_sol(sol)
    n = n0+ν*ξ
    return n
end

function equilibrium_conditions(model::ReactiveModel,F,p,T,n0,ξ,ν,Keq)
    n = n0+ν*ξ

    a = activity(model.eosmodel,p,T,n)

    F[1:end] = sum(log.(a).*ν,dims=1) + log.(Keq)
    return F
end

function x0_equilibrium_conditions(model::ReactiveModel,p,T,n0)
    
    z0 = zeros(length(model.components))
    z0[1] = log10(p)
    return z0
end

export equilibrate