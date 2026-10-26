abstract type LJScalingModel <: EoSModel end

function σϵ_m_vdw1f(ϵ,σ,m,z)
    σϵ_mix = zero(Base.promote_eltype(ϵ,m,σ,z))
    σ_mix = zero(Base.promote_eltype(m,σ,z))
    m̄ = zero(Base.promote_eltype(m,z))
    ∑z = zero(eltype(z))
    for i ∈ eachindex(z)
        _zi = z[i]
        zi,σ3 = _zi*m[i],σ[i,i]^3
        zii = zi*zi
        σϵ_mix += zii*σ3*ϵ[i,i]
        σ_mix += zii*σ3
        m̄ += zi
        ∑z += _zi
        for j ∈ 1:(i-1)
            σ3ij,zij = σ[i,j]^3,zi*z[j]*m[j]
            σϵ_mix += 2*zij*σ3ij*ϵ[i,j]
            σ_mix += 2*zij*σ3ij
        end
    end
    return σ_mix/(m̄*m̄), σϵ_mix/σ_mix, m̄/∑z
end
σϵ_m_vdw1f(ϵ,σ,z) = σϵ_vdw1f(ϵ,σ,z)

function σϵ_vdw1f(ϵ,σ,z)
    σϵ_mix = zero(Base.promote_eltype(ϵ,σ,z))
    σ_mix = zero(Base.promote_eltype(σ,z))
    ∑z = zero(eltype(z))
    for i ∈ eachindex(z)
        zi,σ3 = z[i],σ[i,i]^3
        zii = zi*zi
        σϵ_mix += zii*σ3*ϵ[i,i]
        σ_mix += zii*σ3
        ∑z += zi
        for j ∈ 1:(i-1)
            σ3ij,zij = σ[i,j]^3,zi*z[j]
            σϵ_mix += 2*zij*σ3ij*ϵ[i,j]
            σ_mix += 2*zij*σ3ij
        end
    end
    return σ_mix/(∑z*∑z), σϵ_mix/σ_mix, one(∑z)
end

function σϵ_m_vdw1f(model::EoSModel,V,T,z)
    ϵ = model.params.epsilon.values
    σ = model.params.sigma.values
    m = model.params.segment.values
    return σϵ_m_vdw1f(ϵ,σ,m,z)
end

σϵ_m_vdw1f(ϵ,σ,V,T,z) = σϵ_m_vdw1f(ϵ,σ,z)

struct LJScalingParam{T} <: ParametricEoSParam{T}
    Mw::SingleParam{T}
    sigma::PairParam{T}
    epsilon::PairParam{T}
end

function LJScalingParam(Mw,sigma,epsilon)
    return build_parametric_param(LJScalingParam,Mw,sigma,epsilon)
end

function transform_params_lj(params)
    sigma = params["sigma"]
    sigma.values .*= 1E-10
    return saft_lorentz_berthelot(params)
end

@newmodelsimple LJScaling LJScalingModel LJScalingParam{T}

function LJScaling(params::LJScalingParam{T}) where T
    comps = params.Mw.components
    return LJScaling(comps,params,comps) #use already allocated vector.
end

LJScaling(param::EoSParam) = LJScaling(LJScalingParam(param.Mw,param.sigma,param.epsilon))

function shape_factors(scaling::LJScalingModel,ref,T,z)
    ϵ, σ = scaling.params.epsilon.values, scaling.params.sigma.values
    σ3, ϵ̄, _ = σϵ_m_vdw1f(ϵ,σ,1.0,T,z)
    f = ϵ̄
    h = N_A*σ3
    return f,h
end

inverse_shape_factors(scaling::LJScalingModel,ref,T,z) = shape_factors(scaling,ref,T,z)

function σϵ_m_vdw1f(model::LJScaling,V,T,z)
    ϵ = model.params.epsilon.values
    σ = model.params.sigma.values
    return σϵ_m_vdw1f(ϵ,σ,V,T,z)
end

function transform_params(::Type{LJScaling},params,components)
    sigma = params["sigma"]
    sigma.values .*= 1E-10
    nc = length(components)
    Mw = get!(params,"Mw") do
        SingleParam("Mw",components,zeros(nc),fill(true,nc))
    end
    return saft_lorentz_berthelot(params)
end

"""
    LJScaling <: LJScalingModel
    LJScaling(components;
    userlocations = String[],
    verbose = false)

## Input parameters

- `sigma`: Single Parameter (`Float64`) - particle size `[Å]`
- `epsilon`: Single Parameter (`Float64`) - dispersion energy `[K]`
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g·mol⁻¹]`
- `k`: Pair Parameter (`Float64`) (optional) - Binary Interaction Parameter (no units)
- `l`: Pair Parameter (`Float64`) (optional) - Binary Interaction Parameter (no units)

## Model Parameters

- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g·mol⁻¹]`
- `sigma`: Pair Parameter (`Float64`) - particle size `[m]`
- `epsilon`: Pair Parameter (`Float64`) - dispersion energy `[K]`

## Description

Lennard-Jones corresponding-states scaling model, used together with a Lennard-Jones reference equation of state.

Their shape factors are defined as:
```
σ³ = Σᵢ Σⱼ xᵢ xⱼ σᵢⱼ³
ϵmix = (1/σ³) Σᵢ Σⱼ xᵢ xⱼ σᵢⱼ³ ϵᵢⱼ
σᵢⱼ = (1 - lᵢⱼ)(σᵢ + σⱼ)/2        (Lorentz)
ϵᵢⱼ = (1 - kᵢⱼ)√(ϵᵢ ϵⱼ)           (Berthelot)
f = ϵmix
h = N_A σ³
```
"""
LJScaling

export LJScaling