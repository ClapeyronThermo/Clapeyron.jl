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
        charge = nothing,
        RSPmodel = ConstRSP,
        userlocations = String[],
        RSPmodel_userlocations = String[],
        verbose = false)

## Input models
- `RSPmodel`: Relative Static Permittivity Model

## Description
This function is used to create a Mean Spherical Approximation model. The MSA term gives the excess Helmholtz energy to account for the electrostatic interactions between ions in solution.

The screening length `Γ` is obtained by minimizing the (implicit) Helmholtz energy `Q(Γ)` with a safeguarded Newton method, see [`Clapeyron.msa_obj`](@ref). Derivatives with respect to volume, temperature and composition are obtained with implicit differentiation, so they do not require differentiating through the iterations.

## References
1. Blum, L. (1974). Solution of a model for the solvent-electrolyte interactions in the mean spherical approximation. The Journal of Chemical Physics, 61(5), 2129–2133. [doi:10.1063/1.1682224](https://doi.org/10.1063/1.1682224)
"""
function MSA(solvents,ions; charge = nothing, RSPmodel=ConstRSP, userlocations=String[], RSPmodel_userlocations=String[], verbose=false)
    solvents = format_components(solvents)
    ions = format_components(ions)
    components = vcat(solvents, ions)

    userlocations = normalize_userlocations(userlocations)
    RSPmodel_userlocations = normalize_userlocations(RSPmodel_userlocations)

    references = default_references(MSA)

    init_RSPmodel = @initmodel RSPmodel(solvents,ions,userlocations = RSPmodel_userlocations, verbose = verbose)

    model = MSA(components, init_RSPmodel, references)
    return model
end

default_references(::Type{MSA}) = ["10.1063/1.1682224"]

function a_res(model::MSAModel, V, T, z, iondata)
    return @f(a_MSA,iondata)
end

a_MSA(ionmodel::MSAModel, V, T, z, iondata) = a_MSA(V, T, z, iondata)


function a_MSA(V::Number, T, z, iondata)
    Z, σ, ϵ_r = iondata
    _0 = zero(Base.promote_eltype(V, T, z, Z, σ, ϵ_r))
    if all(iszero,Z)
        return _0
    end
    Γ = screening_length(V, T, z, Z, σ, ϵ_r)
    #Σ1 = ∑zᵢσᵢ³ (independent of Γ)
    #Σ2 = ∑zᵢσᵢZᵢ/(1+Γσᵢ)
    #Σ3 = ∑zᵢσᵢ³/(1+Γσᵢ)
    #Σ4 = ∑zᵢZᵢ²/(1+Γσᵢ)
    ∑z = sum(z)
    ρ = N_A*∑z/V
    Σ1,Σ2,Σ3,Σ4 = msa_sums(Γ, z, Z, σ)
    Δ = 1 - π*ρ/6*Σ1/∑z
    Ω = 1 + π*ρ/(2*Δ)*Σ3/∑z
    Pn = ρ/Ω*Σ2/∑z
    U_MSA = -e_c^2*V/(4π*ϵ_0*ϵ_r)*(Γ*ρ*Σ4/∑z + π/(2Δ)*Ω*Pn^2)
    return (U_MSA + Γ^3*k_B*T*V/(3π))/(N_A*k_B*T*∑z)
end

function screening_length(model::MSAModel,V,T,z,iondata)
    Z, σ, ϵ_r = iondata
    return screening_length(V, T, z, Z, σ, ϵ_r)
end

"""
    screening_length(V, T, z, Z, σ, ϵ_r)

Returns the MSA screening length `Γ` `[1/m]`, as the minimizer of [`Clapeyron.msa_obj`](@ref).

The solution is obtained in the following steps:
1. Initial estimate from the closed-form solution of the MSA closure with `Nᵢ = 0` (no `ησᵢ²` term), using a charge-weighted average ion diameter.
2. Three damped successive substitution steps on the MSA closure ([`Clapeyron.msa_Γ_ss`](@ref)).
3. Upper bound of `Γ` from `(2Γ)² ≤ κ²`.
4. Newton iterations on `Q(Γ)` ([`Clapeyron.msa_newton`](@ref)).
"""
function screening_length(V, T, z, Z, σ, ϵ_r)
    _0 = zero(Base.promote_eltype(V, T, z, Z, σ, ϵ_r))
    data = (V, T, z, Z, σ, ϵ_r)
    λdata = primalval(data)
    λΓ = msa_solver(λdata)
    iszero(λΓ) && return _0
    any(has_dual,data) || return λΓ #no AD: nothing to reconstruct
    return __gradients_for_root_finders(λΓ, data, λdata, msa_ad_cond)
end

#=
Γ-dependent sums over the ions (Zᵢ ≠ 0). `z` are mole amounts, so the number density of ion `i` is `ρᵢ = N_A*zᵢ/V`.
Σ1 = ∑zᵢσᵢ³ (independent of Γ)
Σ2 = ∑zᵢσᵢZᵢ/(1+Γσᵢ)
Σ3 = ∑zᵢσᵢ³/(1+Γσᵢ)
Σ4 = ∑zᵢZᵢ²/(1+Γσᵢ)
=#
function msa_sums(Γ, z, Z, σ)
    _0 = zero(Base.promote_eltype(Γ, z, Z, σ))
    Σ1,Σ2,Σ3,Σ4 = _0,_0,_0,_0
    iions = iions_iterator(Z)
    for i in iions
        zi,Zi,σi = z[i],Z[i],σ[i]
        σi3 = σi*σi*σi
        Σ1 += zi*σi3
        Γσp1 = 1 + Γ*σi
        Σ2 += zi*σi*Zi/Γσp1
        Σ3 += zi*σi3/Γσp1
        Σ4 += zi*Zi*Zi/Γσp1
    end
    return Σ1,Σ2,Σ3,Σ4
end

#η = π/(2Δ)*Pn, with Δ = 1 - π/6*∑ρᵢσᵢ³, Ω = 1 + π/(2Δ)*∑ρᵢσᵢ³/(1+Γσᵢ) and Pn = (1/Ω)*∑ρᵢσᵢZᵢ/(1+Γσᵢ)
function msa_η(ρ, Σ1, Σ2, Σ3)
    Δ = 1 - π*ρ/6*Σ1
    Ω = 1 + π*ρ/(2Δ)*Σ3
    Pn = ρ*Σ2/Ω
    return π*Pn/(2Δ)
end

"""
    msa_obj(Γ, data)

Objective function of the MSA screening length minimization:

```julia
Q = V*Γ^3/(3π) - V*e^2/(4π*ϵ_r*ϵ_0*k_B*T) * ∑ρᵢZᵢ*(Γ*Zᵢ + η*σᵢ)/(1 + σᵢΓ)
```

where `ρᵢ = N_A*zᵢ/V`, and `η` is a function of `Γ`.
`Q` is dimensionless, and equal to `A_MSA/(k_B*T)` at the solution, that is, the reduced MSA Helmholtz energy is `Q(Γ)/(N_A*∑z)`.
"""
function msa_obj(Γ, data)
    V, T, z, Z, σ, ϵ_r = data
    ρ = N_A/V
    Σ1,Σ2,Σ3,Σ4 = msa_sums(Γ, z, Z, σ)
    η = msa_η(ρ, Σ1, Σ2, Σ3)
    c = e_c^2/(4π*ϵ_0*ϵ_r*k_B*T)
    return V*Γ*Γ*Γ/(3π) - V*c*ρ*(Γ*Σ4 + η*Σ2)
end

#solves ∂Q/∂Γ = 0, for primal values of `data`.
function msa_solver(data; kwargs...)
    V, T, z, Z, σ, ϵ_r = data
    κ = debye_length(V, T, z, ϵ_r, Z)
    iszero(κ) && return zero(κ)
    Γmax = κ/2 #(2Γ)² ≤ κ²
    #initial estimate:
    #Nᵢ = 0 gives Γ² = π*c*∑ρᵢZᵢ²/(1+Γσᵢ)². With a single effective diameter, this is Γ*(1 + Γσ) = κ/2.
    Σw,Σwσ = zero(κ),zero(κ)
    iions = iions_iterator(Z)
    for i in iions
        w = z[i]*Z[i]*Z[i]
        Σw += w
        Σwσ += w*σ[i]
    end
    σeff = Σwσ/Σw
    Γ = κ/(1 + sqrt(1 + 2*κ*σeff))
    #damped successive substitution
    for _ in 1:3
        Γ = (Γ + msa_Γ_ss(Γ, data))/2
    end
    Γ = min(Γ, Γmax)
    f = Base.Fix2(msa_obj, data)
    return msa_newton(f, Γ, Γmax; kwargs...)
end

"""
    msa_Γ_ss(Γ, data)

One successive substitution step of the MSA closure: `Γ² = π*e²/(4π*ϵ_r*ϵ_0*k_B*T) * ∑ρᵢ*((Zᵢ - ησᵢ²)/(1 + Γσᵢ))²`.
Returns the right-hand side as the new value of `Γ`.
"""
function msa_Γ_ss(Γ, data)
    V, T, z, Z, σ, ϵ_r = data
    ρ = N_A/V
    Σ1,Σ2,Σ3,_ = msa_sums(Γ, z, Z, σ)
    η = msa_η(ρ, Σ1, Σ2, Σ3)
    Σ = zero(η)
    iions = iions_iterator(Z)
    for i in iions
        σi = σ[i]
        Ni = (Z[i] - η*σi*σi)/(1 + Γ*σi)
        Σ += z[i]*Ni*Ni
    end
    c = e_c^2/(4π*ϵ_0*ϵ_r*k_B*T)
    return sqrt(π*c*ρ*Σ)
end

#implicit AD condition.
msa_ad_cond(Γ, data) = msa_Γ_ss(Γ, data) - Γ

"""
    msa_newton(f, Γ0, Γmax; rtol = eps(one(Γ0))^(4/5), atol = zero(Γ0), max_iters = 100)

Safeguarded Newton minimization of the one-dimensional function `f(Γ)` in `[0,Γmax]`.

"""
function msa_newton(f::F, Γ0::T, Γmax::T; rtol = eps(one(T))^(4/5), atol = zero(T), max_iters = 100) where {F,T<:Real}
    
    Γ = Γ0
    Γlo = zero(T)
    Γhi = Γmax
    for it in 1:max_iters
        _, ∂Q, ∂²Q = Solvers.f∂f∂2f(f, Γ)
        isnan(Γ) && return Γ
        Γnew = Γ - ∂Q/∂²Q
        abs(Γnew - Γ) <= max(atol, rtol*abs(Γnew)) && return Γnew
        if ∂Q > 0 #Γ is over the minimizer
            Γhi = Γ
        else
            Γlo = Γ
        end
        #also catches NaNs.
        if !(Γlo < Γnew < Γhi)
            Γnew = (Γlo + Γhi)/2
        end
        Γ = Γnew
    end
    return T(NaN)
end