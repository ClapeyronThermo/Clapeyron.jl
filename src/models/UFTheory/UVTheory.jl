abstract type UVModel <: EoSModel end
abstract type UVWCAModel <: UVModel end

struct UVParam <: EoSParam
    Mw::SingleParam{Float64}
    sigma::SingleParam{Float64}
    epsilon::SingleParam{Float64}
end


is_splittable(::UVModel) = false

@newmodel UVWCA UVWCAModel UVParam

UVWCA

export UVWCA
function UVWCA(components;
    idealmodel=BasicIdeal,
    userlocations=String[],
    ideal_userlocations=String[],
    verbose=false,)
    params = getparams(components, ["UVUF/UV","properties/molarmass.csv"]; userlocations=userlocations, verbose=verbose)
    
    Mw = params["Mw"]
    sigma = params["sigma"]
    epsilon = params["epsilon"]

    packagedparams = UVParam(Mw, sigma, epsilon)
    references = ["10.1021/ie0003887", "10.1021/ie010954d"]

    model = UVWCA(packagedparams, idealmodel; ideal_userlocations, references, verbose)
    return model
end

function a_res(model::UVWCAModel, V, T, z)
    _data = @f(data)
    (ρ̄,β̄,d̄) = _data
    c = UVWCAconsts.corr_ϕ
    ϕᵤ = tanh(c[1]*ρ̄+c[2]*ρ̄^c[3]+c[4]*ρ̄^c[5])
    return @f(aₕₛ,_data)+@f(Δaᵘ,_data)+(1-ϕᵤ)*(@f(ΔB₂,_data)-@f(ΔBᵘ,_data))*ρ̄
end

function data(model::UVWCAModel,V,T,z)
    σ = model.params.sigma.values[1]
    ϵ = model.params.epsilon.values[1]
    ρ̄  = dot(z,σ.^3)/V
    β̄ = ϵ/T
    d̄ = (2/(1+√(β̄^-1)))^(1/6)
    return (ρ̄,β̄,d̄)
end

function aₕₛ(model::UVWCAModel, V, T, z,_data=@f(data))
    (ρ̄,β̄,d̄) = _data
    η = π/6*ρ̄*d̄^3
    return (4η-3η^2)/(1-η)^2
end

function Δaᵘ(model::UVWCAModel, V, T, z,_data=@f(data))
    (ρ̄,β̄,d̄) = _data
    c = UVWCAconsts.corr_Δaᵘ
    return 2π*ρ̄*β̄*((c[1]+c[2]*d̄^3)+(c[3]+c[4]*d̄^3)*ρ̄+(c[5]+c[6]*d̄^3)*ρ̄^2+(c[7]+c[8]*d̄^3)*ρ̄^3)
end

function ΔBᵘ(model::UVWCAModel, V, T, z,_data=@f(data))
    (ρ̄,β̄,d̄) = _data
    c = UVWCAconsts.corr_Δaᵘ
    return 2π*β̄*(c[1] + c[2]*d̄^3)
end

function ΔB₂(model::UVWCAModel, V, T, z,_data=@f(data))
    (ρ̄,β̄,d̄) = _data
    c = UVWCAconsts.corr_ΔB₂
    T̄ = β̄^-1
    C₁ = c[1] + 4*c[2]
    C₂ = c[3] + 4*c[4] + 16*c[5] + 64*c[6]
    C₃ = 4*c[7] + 16*c[8] + 64*c[9]
    B₂ = 1.345-1.336*β̄-3.85*β̄^2+1.295*β̄^3-0.416*β̄^4
    B₂₀ = 2/3*π*2^(1/2)*(1+√(4π*T̄)+ C₁*T̄ + C₂*T̄^(3/2)+C₃*T̄^2)^(-1/8)
    return B₂-B₂₀
end

const UVWCAconsts = (
    corr_ϕ   = (0.11072527,2.5809656,1.2333385,2.3534453,4.1241290),
    corr_Δaᵘ = (-0.8513603,0.2977649,-0.1287404,-0.3377645,0.7343770,-0.9839062,-0.4592915,0.7963438),
    corr_ΔB₂ = (1.92840364363978e0,4.43165896265079e-1,5.20120816141761e-1,1.82526759234412e-1,1.10319989659929e-2,-7.97813995328348e-5,1.29885156087242e-2,6.41039871789327e-3,1.85866741090323e-5)
)

function lb_volume(model::UVWCAModel,z=SA[1.])
    σ = model.params.sigma.values[1]
    return π/6*σ^3
end

function x0_volume_liquid(model::UVWCAModel,T,z=SA[1.])
    v_lb = lb_volume(model,z)
    return v_lb*1.5
end

function x0_crit_pure(model::UVWCAModel)
    lb_v = lb_volume(model)
    (2.0, log10(lb_v/0.4))
end

function T_scale(model::UVWCAModel,z=SA[1.0])
    ϵ = model.params.epsilon.values[1]
    return ϵ
end

function p_scale(model::UVWCAModel,z=SA[1.0])
    ϵ = model.params.epsilon.values[1]
    σ = model.params.sigma.values[1]
    return ϵ/σ^3
end