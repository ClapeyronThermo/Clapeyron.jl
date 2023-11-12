struct ADPCSAFTParam <: EoSParam
    Mw::SingleParam{Float64}
    segment::SingleParam{Float64}    ##Note: The model was created for spherica molecules (segment = 1). This parameter is not used in computation since it is unecessary, but it is included to not have any issues.
    r_c::SingleParam{Float64}
    sigma::PairParam{Float64}
    epsilon::PairParam{Float64}
    epsilon_assoc::AssocParam{Float64}
    theta_c::AssocParam{Float64}
    c1::SingleParam{Float64}
    c2::SingleParam{Float64}
    c3::SingleParam{Float64}
end

abstract type ADPCSAFTModel <: SAFTModel end
@newmodel ADPCSAFT ADPCSAFTModel ADPCSAFTParam
default_references(::Type{ADPCSAFT}) = ["10.1002/aic.17342"]
default_locations(::Type{ADPCSAFT}) = ["SAFT/ADPCSAFT","properties/molarmass.csv"]
function transform_params(::Type{ADPCSAFT},params, components)
    if length(components) != 1
        "This model was created only for a single component. It has yet to be extended to mixtures" |> ArgumentError |> throw
    end
    sigma = params["sigma"]
    sigma.values .*= 1E-10
    return saft_lorentz_berthelot(params)
end

"""
    ADPCSAFTModel <: SAFTModel
    ADPCSAFT(components; 
    idealmodel=BasicIdeal,
    userlocations=String[],
    ideal_userlocations=String[],
    verbose=false,
    assoc_options = AssocOptions())
    
## Input parameters
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g/mol]`
- `r_c`: Single Parameter (`Float64`)
- `c1`: Single Parameter (`Float64`)
- `c2`: Single Parameter (`Float64`)
- `c3`: Single Parameter (`Float64`)
- `segment`: Single Parameter (`Float64`) - Number of segments (no units)
- `sigma`: Single Parameter (`Float64`) - Segment Diameter [`A°`]
- `epsilon`: Single Parameter (`Float64`) - Reduced dispersion energy  `[K]`
- `k`: Pair Parameter (`Float64`) (optional) - Binary Interaction Paramater (no units)
- `epsilon_assoc`: Association Parameter (`Float64`) - Reduced association energy `[K]`
- `bondvol`: Association Parameter (`Float64`) - Association Volume `[m^3]`
- `theta_c`: Association Parameter (`Float64`)

## Model Parameters
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g/mol]`
- `r_c`: Single Parameter (`Float64`)
- `c1`: Single Parameter (`Float64`)
- `c2`: Single Parameter (`Float64`)
- `c3`: Single Parameter (`Float64`)
- `segment`: Single Parameter (`Float64`) - Number of segments (no units)
- `sigma`: Pair Parameter (`Float64`) - Mixed segment Diameter `[m]`
- `epsilon`: Pair Parameter (`Float64`) - Mixed reduced dispersion energy`[K]`
- `epsilon_assoc`: Association Parameter (`Float64`) - Reduced association energy `[K]`
- `bondvol`: Association Parameter (`Float64`) - Association Volume
- `theta_c`: Association Parameter (`Float64`)

## Input models
- `idealmodel`: Ideal Model

## Description
modified Perturbed-Chain SAFT (PC-SAFT) with an association dependent hard sphere diameter.Currently only works for water.

!!! warning "numerically unstable"
    Due to its functional form, DAPT is numerically unstable. Please use `big` Floats for most calculations.

## References
1. Marshall, B. D. (2021). A modified perturbed chain‐statistical associating fluid theory equation of state for water which includes an association dependent hard sphere diameter. AIChE Journal. American Institute of Chemical Engineers, 67(10). [doi:10.1002/aic.17342](https://doi.org/10.1002/aic.17342)
"""
ADPCSAFT

export ADPCSAFT

recombine_impl!(model::ADPCSAFTModel) = recombine_saft!(model)

function a_res(model::ADPCSAFT, V, T, z)
    _data = @f(data)
    return @f(a_hs,_data) + @f(a_disp,_data) + @f(a_assoc,_data)
end

function data(model::ADPCSAFT, V, T, z)
    Xassoc = @f(X)
    η_a = @f(η_as,Xassoc)
    return Xassoc,η_a
end

function η_as(model::ADPCSAFT, V, T, z, X_ = @f(X))
    σ = model.params.sigma.values[1,1]
    c_1 = model.params.c1.values[1]
    c_2 = model.params.c2.values[1]
    c_3 = model.params.c3.values[1]
    ρ = N_A*∑(z)/V
    η = (π/6)*σ^3*ρ
    XA = X_[1][1]
    X4 = (1.0 - XA)^4
    ψ = 1.0 + c_1*X4 + c_2*X4^2 + c_3*X4^3
    return η*ψ^3
end

function a_hs(model::ADPCSAFTModel, V, T, z,_data = @f(data))
    _,η_a = _data
    ∑z = ∑(z)
    return (4*η_a - 3*η_a^2)/(1-η_a)^2
end

function a_disp(model::ADPCSAFTModel, V, T, z,_data = @f(data))
    _,η_a = _data
    σ = model.params.sigma.values[1,1]
    ρ = N_A*∑(z)/V
    ε = model.params.epsilon.values[1,1]
    I_1 = I_hc(model::ADPCSAFTModel, V, T, z, 1, η_a)
    I_2 = I_hc(model::ADPCSAFTModel, V, T, z, 2, η_a)
    C_1 = (1 + (8*η_a - 2*η_a^2)/(1 - η_a)^4)^-1
    return -2*π*ρ*I_1*σ^3*ε/T - π*ρ*I_2*C_1*σ^3*(ε/T)^2
end

function I_hc(model::ADPCSAFTModel, V, T, z, n,η_a = @f(η_as))
    if n == 1
        corr = ADPCSAFTconsts.corr1
    elseif n == 2
        corr = ADPCSAFTconsts.corr2
    end
    sum = zero(V*T)
    for i = 1:7
        corr1, _ , _ = corr[i] #segment == 1
        ii = i - 1
        sum = sum + corr1*η_a^ii
    end
    return sum
end

function I(model::ADPCSAFTModel, V, T, z, l_c)
    σ = model.params.sigma.values[1,1]
    η = (π/6)*N_A*∑(z)/V*σ^3
    f = evalpoly(η,(3,3,-1))
    k1 = cbrt(2*η*f)
    #a lot of numerical error in this exact expression, rewritten to minimize said error
    k2 = sqrt(evalpoly(η,(9,18,3,-6,3)))/f
    #y_1 = k1*cbrt(k2+1)
    #y_2 = k1*cbrt(k2-1)
    y_1 = k1*EoSFunctions.cbrtp1(k2)
    y_2 = k1*EoSFunctions.cbrtm1(k2)
    z_d = (y_1 - y_2)*k1
    z_s = (y_1 + y_2)*k1
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
    res = (a_1*t_1 + a_2*t_2 + a_3*t_3)
end


function a_assoc(model::ADPCSAFTModel, V, T, z,_data = @f(data))
    X_,_ = _data
    XA = X_[1][1]
    return 4*(log(XA) - XA/2 + 1/2)
end

function X(model::ADPCSAFTModel, V, T, z,_data = nothing)
    σ = model.params.sigma.values[1,1]
    θ_c = model.params.theta_c.values[1,1][2,1]
    κ = (1 - cos(θ_c*π/180))^2/4
    ε_as = model.params.epsilon_assoc.values[1,1][2,1]
    f = exp(ε_as/(T))-1
    ρ = N_A*∑(z)/V
    r_c = model.params.r_c.values[1,1]
    D_AB = 4*π*f*κ*σ^3*@f(I,r_c)
    return @assoc_loop Xold Xnew for i ∈ @comps, a ∈ @sites(i)
        Xnew[i][a] = 1/(1 + 2*ρ*Xold[i][a]*D_AB)
    end
end

const ADPCSAFTconsts = (
    corr1 =
    [(0.9105631445,-0.3084016918, -0.0906148351),
    (0.6361281449, 0.1860531159, 0.4527842806),
    (2.6861347891, -2.5030047259, 0.5962700728),
    (-26.547362491, 21.419793629, -1.7241829131),
    (97.759208784, -65.255885330, -4.1302112531),
    (-159.59154087, 83.318680481, 13.776631870),
    (91.297774084, -33.746922930, -8.6728470368)],

    corr2 =
    [(0.7240946941, -0.5755498075, 0.0976883116),
    (2.2382791861, 0.6995095521, -0.2557574982),
    (-4.0025849485, 3.8925673390, -9.1558561530),
    (-21.003576815, -17.215471648, 20.642075974),
    (26.855641363, 192.67226447, -38.804430052),
    (206.55133841, -161.82646165, 93.626774077),
    (-355.60235612, -165.20769346, -29.666905585)]
)