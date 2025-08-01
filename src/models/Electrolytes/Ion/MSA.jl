abstract type MSAModel <: IonModel end

struct MSA{ϵ} <: MSAModel
    components::Array{String,1}
    RSPmodel::ϵ
    references::Array{String,1}
end

export MSA
"""
    MSA(solvents::Array{String,1},
        ions::Array{String,1};
        RSPmodel = ConstRSP,
        userlocations = String[],
        RSPmodel_userlocations = String[],
        verbose = false)

## Input models
- `RSPmodel`: Relative Static Permittivity Model

## Description
This function is used to create a Mean Spherical Approximation model. The MSA term gives the excess Helmholtz free energy to account for the electrostatic interactions between ions in solution.

## References
1. Blum, L. (1974). Solution of a model for the solvent-electrolyte interactions in the mean spherical approximation. The Journal of Chemical Physics, 61(5), 2129–2133. [doi:10.1063/1.1682224](https://doi.org/10.1063/1.1682224)
"""
function MSA(solvents,ions; RSPmodel=ConstRSP, userlocations=String[], RSPmodel_userlocations=String[], verbose=false)
    components = deepcopy(ions)
    prepend!(components,solvents)

    references = default_references(MSA)

    init_RSPmodel = @initmodel RSPmodel(solvents,ions,userlocations = RSPmodel_userlocations, verbose = verbose)

    model = MSA(components, init_RSPmodel, references)
    return model
end

default_references(::Type{MSA}) = ["10.1063/1.1682224"]

function a_res(model::MSAModel, V, T, z, iondata)
    return @f(a_MSA,iondata)
end

function a_MSA(ionmodel::MSAModel, V, T, z, iondata)
    Z, σ, ϵ_r = iondata
    return a_MSA(V, T, z, Z, σ, ϵ_r)
end

function a_MSA(V::Number, T, z, Z, σ, ϵ_r)
    nc = length(Z)
    iions = @iions
    if all(iszero,Z)
        return zero(Base.promote_eltype(V, T, z, Z, σ, ϵ_r))
    end
    ∑z = sum(z)
    ρ = N_A*sum(z)/V
    Γ = screening_length(V, T, z, Z, σ, ϵ_r)
    Δ = 1-π*ρ/6*sum(z[i]*σ[i]^3 for i ∈ iions)/∑z
    Ω = 1+π*ρ/(2*Δ)*sum(z[i]*σ[i]^3/(1+Γ*σ[i]) for i ∈ iions)/∑z
    Pn = ρ/Ω*sum(z[i]*σ[i]*Z[i]/(1+Γ*σ[i]) for i ∈ iions)/∑z
    U_MSA = -e_c^2*V/(4π*ϵ_0*ϵ_r)*(Γ*ρ*sum(z[i]*Z[i]^2/(1+Γ*σ[i]) for i ∈ iions)/∑z + π/(2Δ)*Ω*Pn^2)
    return (U_MSA+Γ^3*k_B*T*V/(3π))/(N_A*k_B*T*sum(z))
end

function screening_length(model::MSAModel,V,T,z,iondata)
    Z, σ, ϵ_r = iondata
    return screening_length(V, T, z, Z, σ, ϵ_r)
end

function screening_length(V, T, z, Z, σ, ϵ_r)
    iions = @iions
    ρ = N_A/V
    nc = length(Z)
    Δ = 1-π*ρ/6*sum(z[i]*σ[i]^3 for i ∈ 1:nc)
    κ = debye_length(V,T,z,ϵ_r,Z)
    Γold = κ
    _0 = zero(Γold)
    iszero(primalval(Γold)) && return _0
    Γnew = _0
    iter = 1
    tol = oneunit(_0)
    k1 = sqrt(π*e_c^2*ρ/(4π*ϵ_0*ϵ_r*k_B*T))

    #step 1: bounded SS
    while tol>1e-12 && iter < 100
        Ω = 1+π*ρ/(2*Δ)*sum(z[i]*σ[i]^3/(1+Γold*σ[i]) for i ∈ iions)
        Pn = ρ/Ω*sum(z[i]*σ[i]*Z[i]/(1+Γold*σ[i]) for i ∈ iions)
        #Q = @. (Z-σ^2*Pn*(π/(2Δ)))./(1+Γold*σ)
        ∑Q2x = _0
        for i ∈ iions
            Qi = (Z[i]-σ[i]^2*Pn*(π/(2Δ)))/(1+Γold*σ[i])
            ∑Q2x += z[i]*Qi^2
        end
        Γnew = k1*sqrt(∑Q2x)
        tol = abs(1-Γnew/Γold)
        
        Γold = Γnew
    end
    return Γnew
end