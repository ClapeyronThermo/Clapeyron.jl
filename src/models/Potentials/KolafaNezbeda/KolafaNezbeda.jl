struct KolafaNezbedaParam{T} <: EoSParam
    Mw::SingleParam{T}
    sigma::PairParam{T}
    epsilon::PairParam{T}
end

abstract type KolafaNezbedaModel <: EoSModel end
@newmodel KolafaNezbeda KolafaNezbedaModel KolafaNezbedaParam false
default_locations(::Type{KolafaNezbeda}) = ["Potentials/KolafaNezbeda"]
default_references(::Type{KolafaNezbeda}) = ["10.1016/0378-3812(94)02573-G"]

function transform_params(::Type{KolafaNezbeda},params)
    sigma = params["sigma"]
    sigma.values .*= 1E-10
    return saft_lorentz_berthelot(params)
end

"""
    KolafaNezbedaModel <: EoSModel

    KolafaNezbeda(components;
    idealmodel = BasicIdeal,
    userlocations = String[],
    ideal_userlocations = String[],
    reference_state = nothing,
    verbose = false)

## Input parameters
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g·mol⁻¹]`
- `sigma`: Single Parameter (`Float64`) - Lennard-Jones size parameter `[Å]`
- `epsilon`: Single Parameter (`Float64`) - Lennard-Jones energy parameter `[K]`
- `k`: Pair Parameter (`Float64`) (optional) - Binary Interaction Parameter (no units)

## Model Parameters
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g·mol⁻¹]`
- `sigma`: Pair Parameter (`Float64`) - Mixed Lennard-Jones size parameter `[m]`
- `epsilon`: Pair Parameter (`Float64`) - Mixed Lennard-Jones energy parameter `[K]`

## Input models
- `idealmodel`: Ideal Model

## Description

Kolafa-Nezbeda equation of state for the Lennard-Jones fluid.

## References
1. Kolafa, J., & Nezbeda, I. (1994). The Lennard-Jones fluid: An accurate analytic and theoretically-based equation of state. Fluid Phase Equilibria, 100, 1–34. [doi:10.1016/0378-3812(94)02573-G](https://doi.org/10.1016/0378-3812(94)02573-G)
"""
KolafaNezbeda

export KolafaNezbeda

#nonlinear parameter γ of the damping function, eq. (15)
const KN_γ = 1.92907278
const KN_ΔB2 = (0.02459877, 0.0, -7.02181962, 2.90616279, -4.13749995, 0.87361369, 0.43102052, -0.58544978)

#residual term Cᵢⱼ, table 3. one tuple per i, coefficients ordered in j = 2:6
const KN_C0 = (2.01546797, -28.17881636, 28.28313847, -10.42402873)
const KN_C1 = (-19.58371655, 75.62340289, -120.70586598, 93.92740328, -27.37737354)
const KN_C2 = (29.34470520, -112.35356937, 170.64908980, -123.06669187, 34.42288969)
const KN_C4 = (-13.37031968, 65.38059570, -115.09233113, 88.91973082, -25.62099890)

function a_res(model::KolafaNezbedaModel, V, T, z)
    _data = @f(data)
    return @f(a_hs,_data) + @f(a_pert,_data)
end

function data(model::KolafaNezbedaModel,V,T,z)
    ϵ, σ = model.params.epsilon.values, model.params.sigma.values
    σ3,ϵ̄,_ = σϵ_m_vdw1f(ϵ,σ,ones(length(model)),V,T,z)
    T̃ = T/ϵ̄
    ρ̃ = N_A*sum(z)*σ3/V
    d̃ = d_kn(T̃)
    η = π*ρ̃*d̃^3/6
    return (η,ρ̃,T̃,d̃)
end

#hybrid Barker-Henderson hard sphere diameter
function d_kn(T̃)
    sT = sqrt(T̃)
    return 1.080142248 + 0.011117524/T̃ - 0.076383859/sT + 0.000693129*sT - 0.063920968*log(T̃)
end

#residual second virial coefficien
function ΔB2_kn(T̃)
    return evalpoly(1/sqrt(T̃),KN_ΔB2)
end

#hard sphere term, associated to the Boublik-Nezbeda hard sphere EoS
function a_hs(model::KolafaNezbedaModel, V, T, z,_data = @f(data))
    η,ρ̃,T̃,d̃ = _data
    return 5/3*log1p(-η) + η*(34 - 33*η + 4*η*η)/(6*(1-η)^2)
end

function a_pert(model::KolafaNezbedaModel, V, T, z,_data = @f(data))
    return @f(a_B2,_data) + @f(a_corr,_data)
end

function a_B2(model::KolafaNezbedaModel, V, T, z,_data = @f(data))
    η,ρ̃,T̃,d̃ = _data
    return ρ̃*ΔB2_kn(T̃)*exp(-KN_γ*ρ̃*ρ̃)
end

function a_corr(model::KolafaNezbedaModel, V, T, z,_data = @f(data))
    η,ρ̃,T̃,d̃ = _data
    u = 1/sqrt(T̃)
    C0 = evalpoly(ρ̃,KN_C0)
    C1 = evalpoly(ρ̃,KN_C1)
    C2 = evalpoly(ρ̃,KN_C2)
    C4 = evalpoly(ρ̃,KN_C4)
    return ρ̃*ρ̃*u*u*(C0 + u*(C1 + u*(C2 + u*u*C4)))
end

function lb_volume(model::KolafaNezbedaModel,T,z)
    ϵ, σ = model.params.epsilon.values, model.params.sigma.values
    σ3,ϵ̄,_ = σϵ_m_vdw1f(ϵ,σ,ones(length(model)),1.0,1.0,z)
    d̃ = d_kn(T/ϵ̄)
    return sum(z)*N_A*σ3*d̃^3*π/6
end

lb_volume(model::KolafaNezbedaModel,z) = lb_volume(model,T_scale(model,z),z)

function T_scale(model::KolafaNezbedaModel,z)
    ϵ, σ = model.params.epsilon.values, model.params.sigma.values
    σ3,ϵ̄,_ = σϵ_m_vdw1f(ϵ,σ,ones(length(model)),1.0,1.0,z)
    return ϵ̄
end

function p_scale(model::KolafaNezbedaModel,z)
    ϵ, σ = model.params.epsilon.values, model.params.sigma.values
    σ3,ϵ̄,_ = σϵ_m_vdw1f(ϵ,σ,ones(length(model)),1.0,1.0,z)
    return R̄*ϵ̄/(N_A*σ3)
end

function x0_crit_pure(model::KolafaNezbedaModel,z)
    ϵ, σ = model.params.epsilon.values, model.params.sigma.values
    σ3,ϵ̄,_ = σϵ_m_vdw1f(ϵ,σ,ones(length(model)),1.0,1.0,z)
    return (1.32, log10(N_A*σ3/0.31))
end

function set_k!(model::KolafaNezbedaModel,k)
    epsilon = model.params.epsilon
    epsilon = epsilon_LorentzBerthelot!(epsilon,k)
    return model
end
