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
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g/mol]`
- `segment`: Single Parameter (`Float64`) - Number of segments (no units)
- `sigma`: Single Parameter (`Float64`) - Segment Diameter [`A°`]
- `epsilon`: Single Parameter (`Float64`) - Reduced dispersion energy  `[K]`
- `lambda`: Single Parameter (`Float64`) - Soft Well range parameter (no units)
- `k`: Pair Parameter (`Float64`) (optional) - Binary Interaction Paramater (no units)
- `l`: Pair Parameter (`Float64`) (optional) - Binary Interaction Paramater (no units)
- `epsilon_assoc`: Association Parameter (`Float64`) - Reduced association energy `[K]`
- `bondvol`: Association Parameter (`Float64`) - Association Volume `[m^3]`

## Model Parameters
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g/mol]`
- `segment`: Single Parameter (`Float64`) - Number of segments (no units)
- `sigma`: Pair Parameter (`Float64`) - Mixed segment Diameter `[m]`
- `lambda`: Pair Parameter (`Float64`) - Mixed Soft Well range parameter (no units)
- `epsilon`: Pair Parameter (`Float64`) - Mixed reduced dispersion energy`[K]`
- `epsilon_assoc`: Association Parameter (`Float64`) - Reduced association energy `[K]`
- `bondvol`: Association Parameter (`Float64`) - Association Volume

## Input models
- `idealmodel`: Ideal Model

## Description

SAFT, Variable Range (VR) ,Square Well (SW)

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
    return model
end

function get_k(model::SAFTVRSWModel)   
    return get_k_geomean(model.params.epsilon)
end

function get_l(model::SAFTVRSWModel)   
    return get_k_mean(model.params.sigma)
end

function a_res(model::SAFTVRSWModel, V, T, z)
    return @f(a_mono) + @f(a_chain) + @f(a_assoc)
end

function a_mono(model::SAFTVRSWModel, V, T, z)
    return @f(a_hs) + @f(a_disp)
end

function a_disp(model::SAFTVRSWModel, V, T, z)
    return @f(a_1) + @f(a_2)
end

function a_hs(model::SAFTVRSWModel, V, T, z)
    ζ0 = @f(ζn,0)
    ζ1 = @f(ζn,1)
    ζ2 = @f(ζn,2)
    ζ3 = @f(ζn,3)
    Σz = sum(z)
    m = model.params.segment.values
    m̄ = dot(z, m)/Σz
    return m̄*bmcs_hs(ζ0,ζ1,ζ2,ζ3)
end

function ζn(model::SAFTVRSWModel, V, T, z, n)
    σ = model.params.sigma.values
    return π/6*@f(ρ_S)*∑(@f(x_S,i)*σ[i,i]^n for i ∈ @comps)
end

function ρ_S(model::SAFTVRSWModel, V, T, z)
    Σz = sum(z)
    m = model.params.segment.values
    m̄ = dot(z, m)
    N = N_A
    return N/V*m̄
end

function x_S(model::SAFTVRSWModel, V, T, z, i)
    Σz = sum(z)
    m = model.params.segment.values
    m̄ = dot(z, m)
    return z[i]*m[i]/m̄
end

function ζ_X(model::SAFTVRSWModel, V, T, z)
    comps = @comps
    σ = model.params.sigma.values
    return π/6*@f(ρ_S)*∑(@f(x_S,i)*@f(x_S,j)*σ[i,j]^3 for i ∈ comps for j ∈ comps)
end

function a_1(model::SAFTVRSWModel, V, T, z)
    comps = @comps
    Σz = sum(z)
    m = model.params.segment.values
    m̄ = dot(z, m)/Σz
    return -m̄/T*@f(ρ_S)*∑(@f(x_S,i)*@f(x_S,j)*@f(a_1,i,j) for i ∈ comps for j ∈ comps)
end

function a_1(model::SAFTVRSWModel, V, T, z, i, j)
    ϵ = model.params.epsilon.values
    λ = model.params.lambda.values
    σ = model.params.sigma.values
    αVDWij = 2π*ϵ[i,j]*σ[i,j]^3*(λ[i,j]^3-1)/3
    return αVDWij * @f(gHS_0,i,j)
end

function ζeff_X(model::SAFTVRSWModel, V, T, z, λ)
    A = SAFTVRSWconsts.A
    ζ_X_ = @f(ζ_X)
    return A * [1; λ; λ^2] ⋅ [ζ_X_; ζ_X_^2; ζ_X_^3]
end

function gHS_0(model::SAFTVRSWModel,V, T, z, i, j)
    λ = model.params.lambda.values
    ζeff_X_ = @f(ζeff_X,λ[i,j])
    return (1-ζeff_X_/2)/(1-ζeff_X_)^3
end

function a_2(model::SAFTVRSWModel, V, T, z)
    comps = @comps
    Σz = sum(z)
    m = model.params.segment.values
    m̄ = dot(z, m)/Σz
    return m̄/T^2*∑(@f(x_S,i)*@f(x_S,j)*@f(a_2,i,j) for i ∈ comps for j ∈ comps)
end

function a_2(model::SAFTVRSWModel, V, T, z, i, j)
    ϵ = model.params.epsilon.values
    ζ0 = @f(ζn,0)
    ζ1 = @f(ζn,1)
    ζ2 = @f(ζn,2)
    ζ3 = @f(ζn,3)
    KHS = ζ0*(1-ζ3)^4/(ζ0*(1-ζ3)^2+6*ζ1*ζ2*(1-ζ3)+9*ζ2^3)
    return 1/2*KHS*ϵ[i,j]*@f(∂a_1╱∂ρ_S,i,j)
end

function ∂a_1╱∂ρ_S(model::SAFTVRSWModel, V, T, z, i, j)
    ϵ = model.params.epsilon.values
    λ = model.params.lambda.values
    σ = model.params.sigma.values
    αij = 2π*ϵ[i,j]*σ[i,j]^3*(λ[i,j]^3-1)/3
    ζ_X_ = @f(ζ_X)
    ζeff_X_ = @f(ζeff_X,λ[i,j])
    A = SAFTVRSWconsts.A
    ∂ζeff_X╱∂ζ_X = A * [1; λ[i,j]; λ[i,j]^2] ⋅ [ζ_X_; 2ζ_X_^2; 3ζ_X_^3]
    # ∂ζeff_X╱∂ζ_X = A * [1; λ[i,j]; λ[i,j]^2] ⋅ [1; 2ζ_X_; 3ζ_X_^2]
    return -αij*(@f(ρ_S)*@f(gHS_0,i,j)+(5/2-ζeff_X_)/(1-ζeff_X_)^4*∂ζeff_X╱∂ζ_X)
end

function a_chain(model::SAFTVRSWModel, V, T, z)
    Σz = sum(z)
    m = model.params.segment.values
    return -∑(z[i]*(log(@f(γSW,i))*(m[i]-1)) for i ∈ @comps)/Σz
end

function γSW(model::SAFTVRSWModel,V, T, z, i)
    ϵ = model.params.epsilon.values[i,i]
    return @f(gSW,i,i)*exp(-ϵ/T)
end

function gSW(model::SAFTVRSWModel,V, T, z, i, j)
    ϵ = model.params.epsilon.values
    return @f(gHS,i,j)+ϵ[i,j]/T*@f(g_1,i,j)
end

function gHS(model::SAFTVRSWModel,V, T, z, i, j)
    σ = model.params.sigma.values
    ζ3 = @f(ζn,3)
    D = σ[i]*σ[j]/(σ[i]+σ[j])*∑(@f(x_S,k)*σ[k]^2 for k ∈ @comps)/∑(@f(x_S,k)*σ[k]^3 for k ∈ @comps)
    return 1/(1-ζ3)+3*D*ζ3/(1-ζ3)^2+2*(D*ζ3)^2/(1-ζ3)^3
end

function g_1(model::SAFTVRSWModel,V, T, z, i, j)
    λ = model.params.lambda.values
    ζ_X_ = @f(ζ_X)
    ζeff_X_ = @f(ζeff_X,λ[i,j])
    A = SAFTVRSWconsts.A
    ∂ζeff_X╱∂ζ_X = A * [1; λ[i,j]; λ[i,j]^2] ⋅ [1; 2ζ_X_; 3ζ_X_^2]
    ∂ζeff_X╱∂λ = A * [0; 1; 2λ[i,j]] ⋅ [ζ_X_; ζ_X_^2; ζ_X_^3]
    return @f(gHS_0,i,j)+(λ[i,j]^3-1)*(5/2-ζeff_X_)/(1-ζeff_X_)^4*(λ[i,j]/3*∂ζeff_X╱∂λ-ζ_X_*∂ζeff_X╱∂ζ_X)
end

function Δ(model::SAFTVRSWModel, V, T, z, i, j, a, b)
    ϵ_assoc = model.params.epsilon_assoc.values
    κ = model.params.bondvol.values
    g = @f(gSW,i,j)
    return g*(exp(ϵ_assoc[i,j][a,b]/T)-1)*κ[i,j][a,b]
end

const SAFTVRSWconsts = (
    A = [2.25855   -1.50349  0.249434;
    -0.66927  1.40049   -0.827739;
    10.1576   -15.0427   5.30827],
)
