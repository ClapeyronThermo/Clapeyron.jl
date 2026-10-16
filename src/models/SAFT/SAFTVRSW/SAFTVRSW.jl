struct SAFTVRSWParam <: EoSParam
    Mw::SingleParam{Float64}
    segment::SingleParam{Float64}
    sigma::PairParam{Float64}
    lambda::PairParam{Float64}
    epsilon::PairParam{Float64}
    epsilon_assoc::AssocParam{Float64}
    bondvol::AssocParam{Float64}
end

abstract type SAFTVRSWModel <: SAFTModel end
@newmodel SAFTVRSW SAFTVRSWModel SAFTVRSWParam
default_references(::Type{SAFTVRSW}) = ["10.1063/1.473101"]
default_locations(::Type{SAFTVRSW}) = ["SAFT/SAFTVRSW","properties/molarmass.csv"]
function transform_params(::Type{SAFTVRSW},params)
    k = get(params,"k",nothing)
    l = get(params,"l",nothing)
    params["sigma"].values .*= 1E-10
    sigma = sigma_LorentzBerthelot(params["sigma"], l)
    epsilon = epsilon_LorentzBerthelot(params["epsilon"], k)
    lambda = lambda_squarewell(params["lambda"], sigma)
    params["sigma"] = sigma
    params["epsilon"] = epsilon
    params["lambda"] = lambda
    return params
end

export SAFTVRSW

"""
    SAFTVRSWModel <: SAFTModel

    SAFTVRSW(components;
    idealmodel = BasicIdeal,
    userlocations = String[],
    ideal_userlocations = String[],
    reference_state = nothing,
    verbose = false,
    assoc_options = AssocOptions())

## Input parameters
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g·mol⁻¹]`
- `segment`: Single Parameter (`Float64`) - Number of segments (no units)
- `sigma`: Single Parameter (`Float64`) - Segment Diameter `[Å]`
- `epsilon`: Single Parameter (`Float64`) - Reduced dispersion energy `[K]`
- `lambda`: Single Parameter (`Float64`) - Square Well range parameter (no units)
- `k`: Pair Parameter (`Float64`) (optional) - Binary Interaction Parameter (no units)
- `l`: Pair Parameter (`Float64`) (optional) - Binary Interaction Parameter (no units)
- `epsilon_assoc`: Association Parameter (`Float64`) - Reduced association energy `[K]`
- `bondvol`: Association Parameter (`Float64`) - Association Volume `[m³]`

## Model Parameters
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g·mol⁻¹]`
- `segment`: Single Parameter (`Float64`) - Number of segments (no units)
- `sigma`: Pair Parameter (`Float64`) - Mixed segment Diameter `[m]`
- `lambda`: Pair Parameter (`Float64`) - Mixed Square Well range parameter (no units)
- `epsilon`: Pair Parameter (`Float64`) - Mixed reduced dispersion energy`[K]`
- `epsilon_assoc`: Association Parameter (`Float64`) - Reduced association energy `[K]`
- `bondvol`: Association Parameter (`Float64`) - Association Volume `[m³]`

## Input models
- `idealmodel`: Ideal Model

## Description

SAFT, Variable Range (VR), Square Well (SW)

## References
1. Gil-Villegas, A., Galindo, A., Whitehead, P. J., Mills, S. J., Jackson, G., & Burgess, A. N. (1997). Statistical associating fluid theory for chain molecules with attractive potentials of variable range. The Journal of chemical physics, 106(10), 4168–4186. [doi:10.1063/1.473101](https://doi.org/10.1063/1.473101)
"""
SAFTVRSW

function recombine_impl!(model::SAFTVRSWModel)
    sigma = model.params.sigma
    epsilon = model.params.epsilon
    lambda = model.params.lambda
    sigma = sigma_LorentzBerthelot!(sigma)
    epsilon = epsilon_LorentzBerthelot!(epsilon)
    lambda_squarewell!(lambda,sigma)
    recombine_assoc!(model)
    return model
end

function get_k(model::SAFTVRSWModel)   
    return get_k_geomean(model.params.epsilon)
end

function get_l(model::SAFTVRSWModel)   
    return get_k_mean(model.params.sigma)
end

function a_res(model::SAFTVRSWModel, V, T, z)
    _data = @f(data)
    return @f(a_mono,_data) + @f(a_chain,_data) + @f(a_assoc,_data)
end

function data(model::SAFTVRSWModel, V, T, z)
    m̄ = dot(z,model.params.segment)
    _D_gHS = @f(D_gHS,m̄)
    σ = model.params.sigma.values
    ζi = @f(ζ0123,diagvalues(σ))
    _ρ_S = @f(ρ_S,m̄)
    ζₓ = @f(ζ_X,m̄,_ρ_S)
    return (_D_gHS,_ρ_S,ζi,ζₓ,m̄)
end

function a_mono(model::SAFTVRSWModel, V, T, z,_data = @f(data))
    return @f(a_hs,_data) + @f(a_disp,_data)
end

function a_hs(model::SAFTVRSWModel, V, T, z,_data = @f(data))
    _D_gHS,_ρ_S,ζi,ζₓ,m̄ = _data
    ζ0,ζ1,ζ2,ζ3 = ζi
    return m̄*bmcs_hs(ζ0,ζ1,ζ2,ζ3)/sum(z)
end

function ρ_S(model::SAFTVRSWModel, V, T, z, m̄ = dot(z,model.params.segment.values))
    Σz = sum(z)
    m = model.params.segment.values
    N = N_A
    return N/V*m̄
end

function x_S(model::SAFTVRSWModel, V, T, z, i)
    Σz = sum(z)
    m = model.params.segment.values
    m̄ = dot(z, m)
    return z[i]*m[i]/m̄
end

function ζ_X(model::SAFTVRSWModel, V, T, z,m̄ = dot(z,model.params.segment.values),_ρ_S = @f(ρ_S,m̄))
    comps = @comps
    σ = model.params.sigma.values
    m = model.params.segment.values
    m̄⁻¹ = 1/m̄
    ζₓ = zero(first(z)+one(eltype(model)))
    for i in @comps
        x_Si = z[i]*m[i]*m̄⁻¹
        ζₓ += x_Si*x_Si*σ[i,i]^3
        for j in 1:(i-1)
            x_Sj = z[j]*m[j]*m̄⁻¹
            ζₓ += 2*x_Si*x_Sj*σ[i,j]^3
        end
    end
    return ζₓ*π/6*_ρ_S
end

function a_disp(model::SAFTVRSWModel, V, T, z,_data = @f(data))
    _D_gHS,_ρ_S,ζi,ζₓ,m̄ = _data
    ζ0,ζ1,ζ2,ζ3 = ζi
    comps = @comps
    Σz = sum(z)
    m = model.params.segment.values
    ϵ = model.params.epsilon.values
    m̄⁻¹ = 1/m̄
    a₁ = zero(first(z)+one(eltype(model))+T+V)
    a₂ = zero(first(z)+one(eltype(model))+T+V)
    KHS = ζ0*(1-ζ3)^4/(ζ0*(1-ζ3)^2+6*ζ1*ζ2*(1-ζ3)+9*ζ2^3)
    for i in @comps
        x_Si = z[i]*m[i]*m̄⁻¹
        a₁ += x_Si*x_Si*@f(a_1,i,i,ζₓ)
        ϵii = ϵ[i,i]
        a₂ += 0.5*KHS*ϵii*ϵii*x_Si*x_Si*@f(∂a_1╱∂ρ_S╱ϵ,i,i,ζₓ)
        for j in 1:(i-1)
            x_Sj = z[j]*m[j]*m̄⁻¹
            a₁ += 2*x_Si*x_Sj*@f(a_1,i,j,ζₓ)
            ϵij = ϵ[i,j]
            a₂ += KHS*ϵij*ϵij*x_Si*x_Sj*@f(∂a_1╱∂ρ_S╱ϵ,i,j,ζₓ)
        end
    end
    return (-m̄/T*_ρ_S*a₁ + m̄/T^2*a₂)/Σz
end

function a_1(model::SAFTVRSWModel, V, T, z, i, j,ζₓ = @f(ζ_X))
    ϵ = model.params.epsilon.values[i,j]
    λ = model.params.lambda.values[i,j]
    σ = model.params.sigma.values[i,j]
    αVDWij = 2π*ϵ*σ^3*(λ*λ*λ-1)/3
    return αVDWij * @f(gHS_0,λ,ζₓ)
end

function ζeff_X(model::SAFTVRSWModel, V, T, z, λ,ζₓ = @f(ζ_X))
    A = SAFTVRSWconsts.A
    return A * SA[1; λ; λ*λ] ⋅ SA[ζₓ; ζₓ*ζₓ; ζₓ*ζₓ*ζₓ]
end

function gHS_0(model::SAFTVRSWModel,V, T, z, λ,ζₓ = @f(ζ_X),ζeff_X_ = @f(ζeff_X,λ,ζₓ))
    return (1-ζeff_X_/2)/(1-ζeff_X_)^3
end


function a_2(model::SAFTVRSWModel, V, T, z, i, j,_data = @f(data))
    _D_gHS,_ρ_S,ζi,ζₓ,m̄ = _data
    ζ0,ζ1,ζ2,ζ3 = ζi
    ϵ = model.params.epsilon.values
    KHS = ζ0*(1-ζ3)^4/(ζ0*(1-ζ3)^2+6*ζ1*ζ2*(1-ζ3)+9*ζ2^3)
    λ = model.params.lambda.values[i,j]
    return 1/2*KHS*ϵ[i,j]^2*@f(∂a_1╱∂ρ_S╱ϵ,i,j,λ,ζₓ,_ρ_S)
end

function ∂a_1╱∂ρ_S╱ϵ(model::SAFTVRSWModel, V, T, z, i, j,ζₓ = @f(ζ_X),λ = model.params.lambda.values[i,j],_ρ_S = @f(ρ_S),ζeff_X_ = @f(ζeff_X,λ))
    σ = model.params.sigma.values
    αij = 2π*σ[i,j]^3*(λ*λ*λ-1)/3
    A = SAFTVRSWconsts.A
    ∂ζeff_X╱∂ζ_X = A * SA[1; λ; λ*λ] ⋅ SA[ζₓ; 2*ζₓ*ζₓ; 3*ζₓ*ζₓ*ζₓ]
    return -αij*(_ρ_S*@f(gHS_0,λ,ζₓ,ζeff_X_)+(5/2-ζeff_X_)/(1-ζeff_X_)^4*∂ζeff_X╱∂ζ_X)
end

function a_chain(model::SAFTVRSWModel, V, T, z,_data = @f(data))
    _D_gHS,_ρ_S,ζi,ζₓ,m̄ = _data
    ζ0,ζ1,ζ2,ζ3 = ζi
    Σz = sum(z)
    m = model.params.segment.values
    res = zero(V+T+Σz+one(eltype(model)))
    for i in @comps
        res -= z[i]*log(@f(γSW,i,_data))*(m[i]-1)
    end
    return res/Σz
end

function γSW(model::SAFTVRSWModel,V, T, z, i,_data = @f(data))
    ϵ = model.params.epsilon.values[i,i]
    return @f(gSW,i,i,_data)*exp(-ϵ/T)
end

function gSW(model::SAFTVRSWModel,V, T, z, i, j,_data = @f(data))
    ϵ = model.params.epsilon.values
    return @f(gHS,i,j,_data)+ϵ[i,j]/T*@f(g_1,i,j,_data)
end

function gHS(model::SAFTVRSWModel,V, T, z, i, j,_data = @f(data))
    _D_gHS,_ρ_S,ζi,ζₓ,m̄ = _data
    ζ0,ζ1,ζ2,ζ3 = ζi
    σ = model.params.sigma.values
    D = σ[i,i]*σ[j,j]/(σ[i,i]+σ[j,j])*_D_gHS
    return 1/(1-ζ3)+3*D*ζ3/(1-ζ3)^2+2*(D*ζ3)^2/(1-ζ3)^3
end

function D_gHS(model::SAFTVRSWModel,V, T, z, m̄ = dot(z,model.params.segment.values))
    m = model.params.segment.values
    σ = model.params.sigma.values
    Σσ²x_Sᵢ = zero(first(z)+one(eltype(model)))
    Σσ³x_Sᵢ = zero(first(z)+one(eltype(model)))
    m̄⁻¹ = 1/m̄
    for i in @comps
        x_Si = m̄⁻¹*z[i]*m[i]
        σi = σ[i,i]
        σ² = σi*σi
        σ³ = σ²*σi
        Σσ²x_Sᵢ += x_Si*σ²
        Σσ³x_Sᵢ += x_Si*σ³
    end
    return Σσ²x_Sᵢ/Σσ³x_Sᵢ
end

function g_1(model::SAFTVRSWModel,V, T, z, i, j,_data = @f(data))
    _D_gHS,_ρ_S,ζi,ζₓ,m̄ = _data
    ζ0,ζ1,ζ2,ζ3 = ζi
    λ = model.params.lambda.values[i,j]
    ζeff_X_ = @f(ζeff_X,λ)
    A = SAFTVRSWconsts.A
    ∂ζeff_X╱∂ζ_X = A * SA[1.0; λ; λ*λ] ⋅ SA[1; 2ζₓ; 3*ζₓ*ζₓ]
    ∂ζeff_X╱∂λ = A * SA[0.0; 1.0; 2λ] ⋅ SA[ζₓ; ζₓ*ζₓ; ζₓ*ζₓ*ζₓ]
    return @f(gHS_0,λ,ζₓ,ζeff_X_)+(λ^3-1)*(5/2-ζeff_X_)/(1-ζeff_X_)^4*(λ/3*∂ζeff_X╱∂λ-ζₓ*∂ζeff_X╱∂ζ_X)
end

function Δ(model::SAFTVRSWModel, V, T, z, i, j, a, b,_data = @f(data))
    ϵ_assoc = model.params.epsilon_assoc.values
    κ = model.params.bondvol.values
    g = @f(gSW,i,j,_data)
    return g*(expm1(ϵ_assoc[i,j][a,b]/T))*κ[i,j][a,b]
end

const SAFTVRSWconsts = (
    A = SA[2.25855   -1.50349  0.249434;
    -0.66927  1.40049   -0.827739;
    10.1576   -15.0427   5.30827],
)
