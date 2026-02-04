
struct LJRefParam <: EoSParam
    epsilon::PairParam{Float64}
    sigma::PairParam{Float64}
    segment::SingleParam{Float64}
    Mw::SingleParam{Float64}
end

struct LJRef{M} <: EmpiricHelmholtzModel
    components::Vector{String}
    params::LJRefParam
    unscaled_lj::M
    references::Vector{String}
end

export LJRef

"""
    LJRef <: EmpiricHelmholtzModel
    LJRef(components;
    userlocations = String[],
    verbose = false)
## Input parameters
- `sigma`: Single Parameter (`Float64`) - particle size `[Å]`
- `epsilon`: Single Parameter (`Float64`) - dispersion energy `[K]`
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g·mol⁻¹]`
- `k`: Pair Parameter (`Float64`) (optional)
- `sigma` mixing coefficient
## Model Parameters
- `sigma`: Pair Parameter (`Float64`) - particle size `[m]`
- `epsilon`: Pair Parameter (`Float64`) - dispersion energy `[K]`
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g·mol⁻¹]`
## Description
Lennard-Jones Reference equation of state. Valid from 0.5 < T/Tc < 7 and pressures up to p/pc = 500.
```
σᵢⱼ = (σᵢ + σⱼ)/2
ϵᵢⱼ = (1-kᵢⱼ)√(ϵⱼϵⱼ)
σ^3 = Σxᵢxⱼσᵢⱼ^3
ϵ = Σxᵢxⱼϵᵢⱼσᵢⱼ^3/σ^3
τᵢ = 1.32ϵᵢ/T
δᵢ = n(Nₐσᵢ^3)/0.31V
a⁰ᵢ(δ,τ) = log(δᵢ) + 1.5log(τᵢ) - 1.515151515τᵢ + 6.262265814
a⁰(δ,τ,z) = ∑xᵢ(a⁰ᵢ + log(xᵢ))
τ = 1.32ϵ/T
δ = n(Nₐσ^3)/0.31V
aʳ(δ,τ)  = aʳ₁+ aʳ₂ + aʳ₃ + aʳ₄
aʳ₁(δ,τ)  = ∑nᵢδ^(dᵢ)τ^(tᵢ), i ∈ 1:6
aʳ₂(δ,τ)  = ∑nᵢexp(-δ^cᵢ)δ^(dᵢ)τ^(tᵢ), i ∈ 7:12
aʳ₃(δ,τ)  = ∑nᵢexp(-ηᵢ(δ - εᵢ)^2 - βᵢ(τ - γᵢ)^2)δ^(dᵢ)τ^(tᵢ), i ∈ 13:23
```
parameters `n`,`t`,`d`,`c`,`η`,`β`,`γ`,`ε` where obtained via fitting.
!!! warning "Multiple component warning"
    The original model was done with only one component in mind. to support multiple components, a VDW 1-fluid mixing rule (shown above) is implemented, but it is not tested.
## References
1. Thol, M., Rutkai, G., Köster, A., Lustig, R., Span, R., & Vrabec, J. (2016). Equation of state for the Lennard-Jones fluid. Journal of physical and chemical reference data, 45(2), 023101. [doi:10.1063/1.4945000](https://doi.org/10.1063/1.4945000)
"""
LJRef

function LJRef(components;
    userlocations = String[],
    verbose = false)

    components = format_components(components)
    params = getparams(components, ["SAFT/PCSAFT"]; userlocations = userlocations, verbose = verbose)
    Mw = params["Mw"]
    params["sigma"].values .*= 1E-10
    k = get(params,"k",nothing)
    sigma = sigma_LorentzBerthelot(params["sigma"])
    epsilon = epsilon_LorentzBerthelot(params["epsilon"], k)
    segment = params["segment"]
    params = LJRefParam(epsilon,sigma,segment,Mw)
    unscaled_lj = TholLJ()
    references = ["10.1063/1.4945000"]
    return LJRef(components,params,unscaled_lj,references)
end


reduced_a_res(model::LJRef,δ,τ) = reduced_a_res(model.unscaled_lj,δ,τ)
#TODO: better relations? EoSRef was done with one fluid in mind.
#this is technically an unsafe extension.
function _v_scale(model::LJRef,z=SA[1.0])
    σ = model.params.sigma.values
    m = model.params.segment.values

    σ_mix = zero(eltype(z))
    comps = length(model)
    ∑z = sum(z)
    m̄ = dot(m,z)/∑z
    m̄inv = 1/m̄
    for i ∈ 1:comps
        zi = z[i]*m[i]*m̄inv
        zii = zi*zi
        σ3 = σ[i,i]^3
        σ_mix += zii*σ3
        for j ∈ 1:(i-1)
            σ3ij = σ[i,j]^3
            zij = zi*z[j]*m[j]*m̄inv
            σ_mix += 2*zij*σ3ij
        end
    end
    return m̄*N_A*σ_mix/(∑z*∑z)
end

function σϵ_m_vdw1f(model,V,T,z)
    ϵ = model.params.epsilon.values
    σ = model.params.sigma.values
    m = model.params.segment.values

    σϵ_mix = zero(eltype(z))
    σ_mix = zero(eltype(z))
    comps = length(model)
    ∑z = sum(z)
    m̄ = dot(m,z)/∑z
    m̄inv = 1/m̄
    for i ∈ 1:comps
        zi = z[i]*m[i]*m̄inv
        zii = zi*zi
        σ3 = σ[i,i]^3
        σϵ_mix += zii*σ3*ϵ[i,i]
        σ_mix += zii*σ3
        for j ∈ 1:(i-1)
            σ3ij = σ[i,j]^3
            zij = zi*z[j]*m[j]*m̄inv
            σϵ_mix += 2*zij*σ3ij*ϵ[i,j]
            σ_mix += 2*zij*σ3ij
        end
    end
    return σ_mix/(∑z*∑z), σϵ_mix/σ_mix, m̄
end

function VT_scale(model,z = SA[1.0])
    σ3, ϵ, m̄ = σϵ_m_vdw1f(model,1.0,1.0,z)
    return m̄*N_A*σ3, ϵ, m̄
end

function T_scale(model::LJRef,z)
    σ, ϵ , m̄ = VT_scale(model,z)
    return ϵ
end

function lb_volume(model::LJRef,T,z)
    val = π/6*_v_scale(model,z)
    return val
end

function eos_impl(model::LJRef,V,T,z)
    Σz = sum(z)
    α0 = a_ideal(model,V,T,z)
    V0,T0,m̄ = VT_scale(model,z)
    τ = 1.32/(T/T0)
    δ = (Σz*V0)/(0.31*V)
    αr = m̄*reduced_a_res(model,δ,τ)
    x1 = R̄*T*Σz*αr
    x2 = R̄*T*α0
    return x1+x2
end


function a_res(model::LJRef,V,T,z = SA[1.0])
    Σz = sum(z)
    ρ = Σz/V
    V0,T0,m̄ = VT_scale(model,z)
    τ = 1.32/(T/T0)
    δ = (ρ*V0)/0.31
    return  m̄*reduced_a_res(model,δ,τ)
end

function eos_res(model::LJRef,V,T,z = SA[1.0])
    Σz = sum(z)
    ρ = Σz/V
    V0,T0,m̄ = VT_scale(model,z)
    τ = 1.32/(T/T0)
    δ = (ρ*V0)/0.31
    αr = m̄*reduced_a_res(model,δ,τ)
    return R̄*T*Σz*αr
end

function x0_sat_pure(model::LJRef,T)
    σ3, ϵ, m̄  = σϵ_m_vdw1f(model,1.0,1.0,SA[1.0])
    Tc = T_scale(model)
    vl0,vv0 = x0_sat_pure(model.unscaled_lj,T/Tc)
    vl0,vv0
    vl = (m̄*N_A*σ3)*vl0
    vv = (m̄*N_A*σ3)*vv0
    return vl,vv
end

function p_scale(model::LJRef,z)
    rhoc = 1/(_v_scale(model,z))
    Tc = T_scale(model,z)
    return R̄*Tc*rhoc
end

#=
LJ ideal model
=#

struct LJRefIdeal{M} <: IdealModel
    components::Vector{String}
    params::LJRefParam
    unscaled_lj::M
    references::Vector{String}
end

LJRefIdeal(lj::LJRef) = LJRefIdeal(lj.components,lj.params,lj.unscaled_lj,lj.references)
LJRef(lj::LJRefIdeal) = LJRef(lj.components,lj.params,lj.unscaled_lj,lj.references)
idealmodel(model::LJRef) = LJRefIdeal(model)

function LJRefIdeal(components;userlocations = String[], verbose = false)
    lj = LJRef(components;userlocations,verbose)
    return LJRefIdeal(lj)
end

function a_ideal(model::Union{LJRef,LJRefIdeal},V,T,z=SA[1.0])
    ∑z = sum(z)
    ϵ = model.params.epsilon.values
    σ = model.params.sigma.values
    m = model.params.segment.values
    lnΣz = log(∑z)
    res = zero(V+T+first(z))
    ρ = ∑z/V
    for i  ∈ @comps
        mᵢ = m[i]
        τᵢ = 1.32/(T/ϵ[i])
        #δᵢ = (mᵢ*N_A*ρ*σ[i]^3)/0.31
        vᵣᵢ = mᵢ*N_A*σ[i]^3 / 0.31
        zᵢ = z[i]
        a₀ᵢ = reduced_a_ideal(model.unscaled_lj,τᵢ)
        res += zᵢ*a₀ᵢ
        res += xlogx(zᵢ,vᵣᵢ)
    end
    res /= ∑z
    res -= log(V)
    return res
end

"""
    LJRefIdeal <: IdealModel
    LJRef(components;
    userlocations = String[],
    verbose = false)

## Input parameters

- `sigma`: Single Parameter (`Float64`) - particle size `[Å]`
- `epsilon`: Single Parameter (`Float64`) - dispersion energy `[K]`
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g·mol⁻¹]`

## Description

Lennard-Jones Reference equation of state. Ideal Part. Valid from 0.5 < T/Tc < 7 and pressures up to p/pc = 500.


```
τᵢ = 1.32ϵᵢ/T
δᵢ = n(Nₐσᵢ^3)/0.31V
a⁰ᵢ(δ,τ) = log(δᵢ) + 1.5log(τᵢ) - 1.515151515τᵢ + 6.262265814
a⁰(δ,τ,z) = ∑xᵢ(a⁰ᵢ + log(xᵢ))

```

`LJRefIdeal` acts as a wrapper of `LJRef` model, you can access it with `LJRef(model::LJRefIdeal)`.

!!! warning "Mutiple component warning"

    The original model was done with only one component in mind. to support multiple components, a VDW 1-fluid mixing rule (shown above) is implemented, but it is not tested.

## References

1. Thol, M., Rutkai, G., Köster, A., Lustig, R., Span, R., & Vrabec, J. (2016). Equation of state for the Lennard-Jones fluid. Journal of physical and chemical reference data, 45(2), 023101. [doi:10.1063/1.4945000](https://doi.org/10.1063/1.4945000)

"""
LJRefIdeal

export LJRefIdeal

