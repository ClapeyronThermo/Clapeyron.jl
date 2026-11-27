"""
    slle_solubility(model::CompositeModel, p, T)

Calculates the phase boundary for solid-liquid-liquid equilibrium of a **ternary** mixture, at a given temperature `T` and pressure `p`.
Returns a matrix containing the composition of the two liquids phases.

Can only function when solid and liquid models are specified within a `CompositeModel` and when the third component is the solute.
"""
function slle_solubility(model::CompositeModel,p,T)
    p = p*one(eltype(model))
    T = T*one(eltype(model))
    if length(model) != 3
        error("SLLE can only be obtained for ternary systems")
    end
    j = findfirst(==([3]),_mapping_split(model))
    j === nothing && throw(ArgumentError("slle_solubility: no solid formed only by the third component"))
    Δμ = sle_Δμ(model,p,T,j)
    μ_ref = sle_μ_ref(model,p,T)
    x0 = x0_slle_solubility(model,p,T,Δμ)

    f!(F,x) = obj_slle_solubility(F,model,p,T,[exp10(x[1]),1-exp10(x[1])-exp10(x[2]),exp10(x[2])],[exp10(x[3]),1-exp10(x[3])-exp10(x[4]),exp10(x[4])],Δμ,μ_ref)
    results = Solvers.nlsolve(f!,x0,LineSearch(Newton(),Backtracking()))
    sol = exp10.(Solvers.x_sol(results))
    x1 = [sol[1],1-sol[1]-sol[2],sol[2]]
    x2 = [sol[3],1-sol[3]-sol[4],sol[4]]
    return (x1,x2)
end

function obj_slle_solubility(F,model,p,T,x1,x2,Δμ,μ_ref)
    nc = length(model)
    lnxγ1 = sle_lnx_γ(model,p,T,x1,μ_ref)
    lnxγ2 = sle_lnx_γ(model,p,T,x2,μ_ref)
    F[1] = lnxγ1[end] - Δμ
    F[2:nc+1] = lnxγ2 .- lnxγ1
    return F
end

function x0_slle_solubility(model,p,T,Δμ)
    x0 = zeros(length(model),2)
    for i in 1:length(model)-1
        z∞ = ones(length(model))*1e-3
        z∞[i] = 1.0
        z∞ ./= sum(z∞)
        γ∞ = activity_coefficient(model.fluid,p,T,z∞)
        z∞[end] = exp(Δμ-log(γ∞[end]))
        z∞ ./= sum(z∞)
        γ∞ = activity_coefficient(model.fluid,p,T,z∞)
        x = zeros(length(model))
        x[i] = 1.0
        x[end] = exp(Δμ-log(γ∞[end]))
        x[findfirst(x.==0)] = 1/γ∞[findfirst(x.==0)]
        x ./= sum(x)

        # x = x0_lle_init(model.fluid,p,T,z0,x)

        x0[:,i] = x
    end

    x0 = [x0[1,1],x0[3,1],x0[1,2],x0[3,2]]

    return log10.(x0)
end
