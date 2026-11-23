#fluid components (sorted) of solid j and their stoichiometric coefficients
function sle_mapping(model::CompositeModel,j)
    idx,ν = _mapping_split(model)[j],_mapping_fractions(model)[j]
    perm = sortperm(idx)
    return idx[perm],ν[perm]
end

function sle_pure_solid_μ(solid,p,T,j)
    solid_j = isone(length(solid)) ? solid : first(index_reduction(solid,(1:length(solid)) .== j))
    return chemical_potential(solid_j,p,T,SA[1.0])[1]
end

#Δμ/RT of solid j: μ of the solid minus Σνₖμₖ of its fluid components in the liquid reference state given by sle_μ_ref
sle_Δμ(model::CompositeModel,p,T,j) = sle_Δμ(model.solid,model,p,T,j)

function sle_Δμ(solid,model,p,T,j)
    idx,ν = sle_mapping(model,j)
    μl = sum(ν .* reference_chemical_potential(model.fluid,p,T,:pure;phase = :l)[idx])
    μs = sle_pure_solid_μ(solid,p,T,j)
    ref = model.solid_reference_state
    length(ref.a0) >= j && (μs += ref.a0[j] + ref.a1[j]*T)
    return (μs - μl)/(Rgas(model.fluid)*T)
end

sle_μ_ref(model::CompositeModel,p,T) = sle_μ_ref(model.solid,model.fluid,p,T)
sle_μ_ref(solid,fluid,p,T) = reference_chemical_potential(fluid,p,T,reference_chemical_potential_type(fluid);phase = :l)

#ln(xᵢγᵢ) of the liquid
sle_lnx_γ(model::CompositeModel,p,T,x,μ_ref) = log.(x .* activity_coefficient(model.fluid,p,T,x;μ_ref,phase = :l))

sle_T_ref(model) = nothing

"""
    sle_solubility(model::CompositeModel, p, T, z; solute)

Calculates the solubility of each component within a solution of the other components, at a given temperature `T`, pressure `p` and composition `z`.
Returns a matrix containing the composition of the SLE phase boundary for each component. If `solute` is specified, returns only the solubility of the specified component.

Can only function when solid and fluid models are specified within a CompositeModel.
"""
function sle_solubility(model::CompositeModel,p,T,z;solute=nothing,x0=nothing)
    solid_components = component_list(model.solid)
    isnothing(solute) && (solute = solid_components)
    nf = length(model.fluid)
    p = p*one(eltype(model))
    T = T*one(eltype(model))
    μ_ref = sle_μ_ref(model,p,T)
    Tm = sle_T_ref(model.solid)
    sol = zeros(length(solute),nf)
    for i in eachindex(solute)
        j = findfirst(==(solute[i]),solid_components)
        Tm !== nothing && T > Tm[j] && error("Temperature above melting point of $(solute[i])")
        idx,νj = sle_mapping(model,j)
        ν = zeros(nf)
        ν[idx] .= νj
        idx_solv = iszero.(ν)
        any(idx_solv) || (idx_solv[first(idx)] = true)
        Δμ = sle_Δμ(model,p,T,j)
        x0i = isnothing(x0) ? x0_sle_solubility(model,p,T,z,idx_solv,ν,Δμ,μ_ref) : x0
        f!(F,x) = obj_sle_solubility(F,model,p,T,z,idx_solv,ν,exp10(x[1]),Δμ,μ_ref)
        results = Solvers.nlsolve(f!,x0i,LineSearch(Newton(),Backtracking()),NEqOptions(),ForwardDiff.Chunk{1}())
        sol[i,:] .= sle_composition(z,idx_solv,ν,exp10(Solvers.x_sol(results)[1]))
    end
    return length(solute) == 1 ? sol[1,:] : sol
end

function sle_composition(z,idx_solv,ν,s)
    x = zeros(Base.promote_eltype(z,s),length(z))
    for k in eachindex(x)
        x[k] = idx_solv[k] ? z[k] : s*ν[k]
    end
    x ./= sum(x)
    return x
end

function obj_sle_solubility(F,model,p,T,z,idx_solv,ν,s,Δμ,μ_ref)
    x = sle_composition(z,idx_solv,ν,s)
    lnxγ = sle_lnx_γ(model,p,T,x,μ_ref)
    F[1] = sum(ν[k]*lnxγ[k] for k in eachindex(ν) if ν[k] > 0) - Δμ
    return F
end

function x0_sle_solubility(model,p,T,z,idx_solv,ν,Δμ,μ_ref)
    z∞ = sle_composition(z,idx_solv,ν,1e-10)
    lnγ∞ = log.(activity_coefficient(model.fluid,p,T,z∞;μ_ref,phase = :l))
    ∑νlnγ∞ = sum(ν[k]*lnγ∞[k] for k in eachindex(ν) if ν[k] > 0)
    isfinite(∑νlnγ∞) || (∑νlnγ∞ = zero(∑νlnγ∞))
    return [min((Δμ - ∑νlnγ∞)/(sum(ν)*log(10)),log10(0.5))]
end

"""
    sle_solubility_T(model::CompositeModel, z, p=1e5; solute=nothing, x0=nothing)

Calculates the solid-liquid coexistence temperature at composition z. Returns the temperature of the most stable solid complex.

Can only function when solid and fluid models are specified within a CompositeModel.
"""
function sle_solubility_T(model::CompositeModel,z,p=1e5;solute=nothing,x0=nothing)
    solid_components = component_list(model.solid)
    isnothing(solute) && (solute = solid_components)
    p = p*one(eltype(model))
    x = z ./ sum(z)
    Tsol = zeros(length(solute))
    for i in eachindex(solute)
        j = findfirst(==(solute[i]),solid_components)
        x0i = isnothing(x0) ? [sle_T_ref(model.solid)[j]] : x0
        f!(F,T) = obj_sle_solubility_T(F,model,p,T[1],x,j)
        results = Solvers.nlsolve(f!,x0i,LineSearch(Newton(),Backtracking()),NEqOptions(),ForwardDiff.Chunk{1}())
        Tsol[i] = Solvers.x_sol(results)[1]
    end
    return maximum(Tsol)
end

function obj_sle_solubility_T(F,model,p,T,x,j)
    idx,ν = sle_mapping(model,j)
    lnxγ = sle_lnx_γ(model,p,T,x,sle_μ_ref(model,p,T))
    F[1] = sum(ν[k]*lnxγ[idx[k]] for k in eachindex(idx)) - sle_Δμ(model,p,T,j)
    return F
end
