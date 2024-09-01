struct DAPTParam <: EoSParam
    Mw::SingleParam{Float64}
    segment::SingleParam{Float64}  ##Note: The model was created for spherica molecules (segment = 1). This parameter is not used in computation since it is unecessary, but it is included to not have any issues.
    r_c::SingleParam{Float64}
    lambda :: SingleParam{Float64}
    sigma::PairParam{Float64}
    epsilon::PairParam{Float64}
    epsilon_assoc::AssocParam{Float64}
    theta_c::AssocParam{Float64}
end

abstract type DAPTModel <: SAFTModel end
@newmodel DAPT DAPTModel DAPTParam
default_references(::Type{DAPT}) = ["10.1016/j.fluid.2019.112252"]
default_locations(::Type{DAPT}) = ["SAFT/DAPT","properties/molarmass.csv"]
function transform_params(::Type{DAPT},params,components)
    if length(components) != 1
        "This model was created only for a single component. It has yet to be extended to mixtures" |> ArgumentError |> throw
    end
    sigma = params["sigma"]
    sigma.values .*= 1E-10
    return saft_lorentz_berthelot(params)
end
"""
    DAPTModel <: SAFTModel
    DAPT(components; 
    idealmodel = BasicIdeal,
    userlocations = String[],
    ideal_userlocations = String[],
    reference_state = nothing,
    verbose = false,
    assoc_options = AssocOptions())

## Input parameters
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g/mol]`
- `r_c`: Single Parameter (`Float64`)
- `lambda`: Single Parameter (`Float64`)
- `segment`: Single Parameter (`Float64`) - Number of segments (no units)
- `sigma`: Single Parameter (`Float64`) - Segment Diameter [`A°`]
- `epsilon`: Single Parameter (`Float64`) - Reduced dispersion energy  `[K]`
- `k`: Pair Parameter (`Float64`) (optional) - Binary Interaction Paramater (no units)
- `epsilon_assoc`: Association Parameter (`Float64`) - Reduced association energy `[K]`
- `theta_c`: Association Parameter (`Float64`)

## Model Parameters
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g/mol]`
- `r_c`: Single Parameter (`Float64`)
- `lambda`: Single Parameter (`Float64`)
- `segment`: Single Parameter (`Float64`) - Number of segments (no units)
- `sigma`: Pair Parameter (`Float64`) - Mixed segment Diameter `[m]`
- `epsilon`: Pair Parameter (`Float64`) - Mixed reduced dispersion energy`[K]`
- `epsilon_assoc`: Association Parameter (`Float64`) - Reduced association energy `[K]`
- `theta_c`: Association Parameter (`Float64`)

## Input models
- `idealmodel`: Ideal Model

## Description
Doubly-Associated Perturbation Theory model. Currently only works for water.

!!! warning "numerically unstable"
    Due to its functional form, DAPT is numerically unstable. Please use `big` Floats for most calculations.

## References
1. Marshall, B. D. (2019). A doubly associated reference perturbation theory for water. Fluid Phase Equilibria, 500(112252), 112252. [doi:10.1016/j.fluid.2019.112252](https://doi.org/10.1016/j.fluid.2019.112252) 
"""
DAPT
#==
function DAPT(components;
    idealmodel = BasicIdeal, userlocations = String[], ideal_userlocations = String[], verbose = false, assoc_options = AssocOptions())
    params,sites = getparams(components, ["SAFT/DAPT","properties/molarmass.csv"]; userlocations = userlocations, verbose = verbose)
    segment = params["m"]
    k = params["k"]  #Note: this is the kij, not the association volume
    Mw = params["Mw"]
    params["sigma"].values .*= 1E-10
    r_c = params["r_c"]   ### Here maybe I will have to add lambda_LorentzBerthelot for the mixing rules
    lambda = params["lambda"]
    sigma = sigma_LorentzBerthelot(params["sigma"])
    epsilon = epsilon_LorentzBerthelot(params["epsilon"], k)
    epsilon_assoc = params["epsilon_assoc"]
    theta_c = params["theta_c"]
    packagedparams = DAPTParam(Mw, segment, r_c, lambda, sigma, epsilon, epsilon_assoc, theta_c)
    references = ["10.1021/ie0003887", "10.1021/ie010954d"]
    model = DAPT(packagedparams, sites, idealmodel; ideal_userlocations = ideal_userlocations, references=references, verbose = verbose)
    return model
end
==#
#DAPT

export DAPT

recombine_impl!(model::DAPTModel) = recombine_saft!(model)

function a_res(model::DAPTModel, V, T, z)
    _data = @f(data)
    return @f(a_hs) + @f(a_disp,_data) + @f(a_assoc,_data)
end

function data(model::DAPTModel,V,T,z)
    Irc = @f(I)
    Xassoc = X(model, V, T, z, Irc)
    return Irc,Xassoc
end

function a_disp(model::DAPTModel, V, T, z,_data = @f(data))
    Irc,X_ = _data
    σ = model.params.sigma.values[1,1]
    ρ = N_A*∑(z)/V
    λ = model.params.lambda.values[1,1]
    ε = model.params.epsilon.values[1,1]
    XA = X_[1][1]
    X4 = (1-XA)^4
    I_ref = @f(I,λ) + X4*(1/(π*ρ*σ^3) - Irc)
    return -2*π*ε*ρ*σ^3*I_ref/(T)
end

function a_hs(model::DAPTModel, V, T, z)
    ∑z = ∑(z)
    σ = model.params.sigma.values[1,1]
    ρ = N_A*∑(z)/V
    η = (π/6)*ρ*σ^3
    return (4*η - 3*η^2)/(1-η)^2
end

function I(model::DAPTModel, V, T, z, l_c = model.params.r_c.values[1])
    ∑z = ∑(z)
    σ = model.params.sigma.values[1,1]
    η = (π/6)*N_A*∑z/V*σ^3
    f = evalpoly(η,(3,3,-1))
    k1 = cbrt(2*η*f)
    #a lot of numerical error in this exact expression, rewritten to minimize said error
    k2 = sqrt(evalpoly(η,(9,18,3,-6,3)))/f

    #y_1 = k1*cbrt(k2+1)
    #y_2 = k1*cbrt(k2-1)
    y_1 = k1*EoSFunctions.cbrtp1(k2)
    y_2 = k1*EoSFunctions.cbrtm1(k2)

    z_d = y_1 - y_2
    z_s = y_1 + y_2
    A = (-2*η + z_d)/(1-η)
    B = (-2*η - 0.5*z_d)/(1-η)
    C = (sqrt(3)*z_s)/(2*(1-η))

    a_1 = (-2*η*(1-η - 3*η^2) + (1 - 3*η - 4*η^2)*z_d + (1+η/2)*z_d^2)/(3*(2*η^2+ z_d^2)*(1-η)^2)
    a_2 = (η*(2+ 4*η - 3*η^2) - (1 - 3*η - 4*η^2)*z_d + 2*(1+η/2)*z_d^2)/(3*(2*η^2+ z_d^2)*(1-η)^2)
    a_3 = ((1-3*η - 4*η^2)*(4*η^2 +z_d^2)+ η*z_d*(2-5*η^2))/(sqrt(3)*z_s*(2*η^2+z_d^2)*(1-η)^2)

    γ = l_c - 1
    t_1 = (exp(A*γ)*(l_c*A-1)-A+1)/A^2
    t_2 = (l_c*exp(B*γ)*(B*cos(C*γ)+C*sin(C*γ))-B)/(B^2+C^2) - (exp(B*γ)*((B^2-C^2)*cos(C*γ) + 2*B*C*sin(C*γ))-(B^2-C^2))/(B^2 + C^2)^2
    t_3 = (l_c*exp(B*γ)*(B*sin(C*γ)-C*cos(C*γ))+C)/(B^2+C^2) - (exp(B*γ)*((B^2-C^2)*sin(C*γ) - 2*B*C*cos(C*γ))+2*B*C)/(B^2 + C^2)^2

    return (a_1*t_1 + a_2*t_2 + a_3*t_3)
end

function a_assoc(model::DAPTModel, V, T, z, _data = @f(data))
    Irc,X_ = _data
    θ_c = model.params.theta_c.values[1,1][2]
    κ = (1 - cos(θ_c*π/180))^2/4
    σ = model.params.sigma.values[1]
    ε = model.params.epsilon_assoc.values[1,1][2]
    f = exp(ε/(T))-1
    ∑z = ∑(z)
    ρ = N_A*∑z/V
    XA = X_[1][1]
    X4 = (1-XA)^4
    Δc0_N = 16*π*f*κ*σ^3*ρ*XA^2*(Irc*(1-X4) + X4/(π*ρ*σ^3))
    return 4*(log(XA) + 1 - XA) - Δc0_N
end

function X(model::DAPTModel, V, T, z,Irc = @f(I))
    σ = model.params.sigma.values[1][1]
    θ_c = model.params.theta_c.values[1,1][2,1]
    κ = (1 - cos(θ_c*π/180))^2/4
    ε_as = model.params.epsilon_assoc.values[1,1][2,1]
    f = exp(ε_as/(T))-1
    ρ = N_A*∑(z)/V
    Xsol = @assoc_loop Xold Xnew for i ∈ @comps, a ∈ @sites(i)
            X4 = (1-Xold[i][a])^4
            c_A = 8*π*κ*σ^3*f*(ρ*Xold[i][a]*(Irc*(1-X4) + X4/(π*ρ*σ^3)) + 2*ρ*(Xold[i][a]^2)*((1 - Xold[i][a])^3)*(Irc - 1/(π*ρ*σ^3)) )
            Xnew[i][a] =1/(1+c_A)
    end
    return Xsol
end

x0_volume_gas(model::DAPTModel,p,T,z) = Rgas(model)*T/p