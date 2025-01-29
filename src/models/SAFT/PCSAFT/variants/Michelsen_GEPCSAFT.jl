
abstract type Michelsen_GEPCSAFTModel <: SAFTModel end

const Michelsen_GEPCSAFTParam = PCSAFTParam

struct Michelsen_GEPCSAFTi{T <: IdealModel,γ} <: Michelsen_GEPCSAFTModel
    components::Array{String,1}
    sites::SiteParam
    activity::γ
    params::Michelsen_GEPCSAFTParam
    idealmodel::T
    assoc_options::AssocOptions
    references::Array{String,1}
end

"""
    Michelsen_GEPCSAFT <: SAFTModel

    Michelsen_GEPCSAFT(components;
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
- `k`: Pair Parameter (`Float64`) (optional) - Binary Interaction Paramater (no units)
- `epsilon_assoc`: Association Parameter (`Float64`) - Reduced association energy `[K]`
- `bondvol`: Association Parameter (`Float64`) - Association Volume `[m^3]`

## Model Parameters
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g/mol]`
- `segment`: Single Parameter (`Float64`) - Number of segments (no units)
- `sigma`: Pair Parameter (`Float64`) - Mixed segment Diameter `[m]`
- `epsilon`: Pair Parameter (`Float64`) - Mixed reduced dispersion energy`[K]`
- `epsilon_assoc`: Association Parameter (`Float64`) - Reduced association energy `[K]`
- `bondvol`: Association Parameter (`Float64`) - Association Volume

## Input models
- `idealmodel`: Ideal Model
- `activity`: Activity model

## Description

Perturbed-Chain SAFT (PC-SAFT), with Gᴱ mixing rule (using the Michelsen limit).

## References
"""
Michelsen_GEPCSAFT

export Michelsen_GEPCSAFT
function Michelsen_GEPCSAFT(components;
    idealmodel = BasicIdeal,
    activity = UNIFAC,
    userlocations = String[],
    ideal_userlocations = String[],
    activity_userlocations = String[],
    assoc_options = AssocOptions(),
    reference_state = nothing,
    verbose = false)

    
    params = getparams(components, ["SAFT/PCSAFT/PCSAFT_like.csv","SAFT/PCSAFT/PCSAFT_unlike.csv","SAFT/PCSAFT/PCSAFT_assoc.csv"]; userlocations = userlocations, verbose = verbose)

    sites = params["sites"]
    segment = params["segment"]
    k = get(params,"k",nothing)
    Mw = params["Mw"]
    params["sigma"].values .*= 1E-10
    sigma = sigma_LorentzBerthelot(params["sigma"])
    epsilon = epsilon_LorentzBerthelot(params["epsilon"], k)
    epsilon_assoc = params["epsilon_assoc"]
    bondvol = params["bondvol"]
    packagedparams = Michelsen_GEPCSAFTParam(Mw, segment, sigma, epsilon, epsilon_assoc, bondvol)

    init_idealmodel = init_model(idealmodel,components,ideal_userlocations,verbose)
    init_activity = init_model(activity,components,activity_userlocations,verbose)
    references = String["10.1021/acs.iecr.2c03464"]
    
    model = Michelsen_GEPCSAFTi(format_components(components),sites,init_activity,packagedparams,init_idealmodel,assoc_options,references)
    set_reference_state!(model,reference_state;verbose)
    return model 
end

function a_res(model::Michelsen_GEPCSAFTModel, V, T, z)
    _data = @f(data)
    return @f(a_hc,_data) + @f(a_disp,_data) + @f(a_assoc,_data)
end

function data(model::Michelsen_GEPCSAFTModel,V,T,z)
    _d = @f(d)
    ζ0,ζ1,ζ2,ζ3 = @f(ζ0123,_d)
    m = model.params.segment.values
    m̄ = dot(z, m)/sum(z)
    return (_d,ζ0,ζ1,ζ2,ζ3,m̄)
end

function a_hc(model::Michelsen_GEPCSAFTModel, V, T, z , _data = @f(data))
    _,_,_,_,η,m̄ = _data
    g_hs = (1-η/2)/(1-η)^3
    a_hs = (4η-3η^2)/(1-η)^2
    return m̄*a_hs - (m̄-1)*log(g_hs)
end

function g_hs(model::Michelsen_GEPCSAFTModel, V, T, z,_data = @f(data))
    _,_,_,_,η,_ = _data
    return (1-η/2)/(1-η)^3
end

function a_hs(model::Michelsen_GEPCSAFTModel, V, T, z)
    _,_,_,_,η,_ = _data
    return (4η-3η^2)/(1-η)^2
end

function Δ(model::Michelsen_GEPCSAFTModel, V, T, z, i, j, a, b,_data = @f(data))
    ϵ_associjab = model.params.epsilon_assoc.values[i,j][a,b]
    κijab = model.params.bondvol.values[i,j][a,b]
    σij = model.params.sigma.values[i,j]
    g_hs_ = @f(g_hs,_data)
    return g_hs_*σij^3*(exp(ϵ_associjab/T)-1)*κijab
end

function a_disp(model::Michelsen_GEPCSAFTModel, V, T, z,_data=@f(data))
    di,ζ0,ζ1,ζ2,ζ3,m̄ = _data
    Σz = sum(z)
    m2ϵσ3₁,m2ϵσ3₂ = @f(m2ϵσ3,_data)
    return -2*π*N_A*Σz/V*@f(I,1,_data)*m2ϵσ3₁ - π*m̄*N_A*Σz/V*@f(C1,_data)*@f(I,2,_data)*m2ϵσ3₂
end

function d(model::Michelsen_GEPCSAFTModel, V, T, z)
    ϵᵢᵢ = diagvalues(model.params.epsilon)
    σᵢᵢ = diagvalues(model.params.sigma)
    return σᵢᵢ .* (1 .- 0.12 .* exp.(-3ϵᵢᵢ ./ T))
end


function C1(model::Michelsen_GEPCSAFTModel, V, T, z,_data=@f(data))
    _,_,_,_,η,m̄ = _data
    return (1 + m̄*(8η-2η^2)/(1-η)^4 + (1-m̄)*(20η-27η^2+12η^3-2η^4)/((1-η)*(2-η))^2)^-1
end

function m2ϵσ3(model::Michelsen_GEPCSAFTModel, V, T, z,_data=@f(data))
    z = [z[i]/sum(z) for i ∈ eachindex(z)]
    d,_,_,_,_,m̄ = _data
    m_i = model.params.segment.values
    m_mix = sum([m_i[i] * z[i] for i ∈ eachindex(z)])
    ϵ_ii = Clapeyron.diagvalues(model.params.epsilon.values)
    σ_ii = Clapeyron.diagvalues(model.params.sigma.values)
    σ_mix = sum([σ_ii[i] * z[i] for i ∈ eachindex(z)])
    gₑ = excess_gibbs_free_energy(model.activity,V,T,z)/(R̄*T)
    bi = [Clapeyron.N_A*m_i[i]*d[i]^3 for i ∈ eachindex(z)]
    b = sum([bi[i]*z[i] for i ∈ eachindex(z)])
    C₂ = gₑ + sum([z[i]*log(b/bi[i]) for i ∈ eachindex(z)])
    α_ii = [m_i[i] * ϵ_ii[i] for i ∈ eachindex(z)]
    α_bar = sum([α_ii[i]*z[i] for i ∈ eachindex(z)])
    a = [-162.1526266, -8.76E-04, -1.49E-10, -13.28867758, -4.46E-01, -1.19E-04, -556.4748244, -471.4372287, 67201.59286]
    Q(α_guess) = a[1] + a[2]*α_guess + a[3]*α_guess^2 + a[4]*log(b) + a[5]*(log(b))^2 + a[6]*log(b)*α_guess + a[7]/α_guess + a[8]/(log(b)) + a[9]/(α_guess^2)
    f(α_guess) = C₂/Q(α_guess) - (α_guess - α_bar)
    α_mix = Clapeyron.Solvers.ad_newton(f, α_bar)
    @show α_mix
    m2ϵσ3₁ = α_mix * m_mix * σ_mix^3
    m2ϵσ3₂ =  α_mix^2*σ_mix^3
    return m2ϵσ3₁,m2ϵσ3₂
end

function I(model::Michelsen_GEPCSAFTModel, V, T, z, n , _data=@f(data))
    _,_,_,_,η,m̄ = _data
    if n == 1
        corr = PCSAFTconsts.corr1
    elseif n == 2
        corr = PCSAFTconsts.corr2
    end
    res = zero(η)
    @inbounds for i ∈ 1:7
        ii = i-1
        corr1,corr2,corr3 = corr[i]
        ki = corr1 + (m̄-1)/m̄*corr2 + (m̄-1)/m̄*(m̄-2)/m̄*corr3
        res +=ki*η^ii
    end
    return res
end

function Ī(model::Michelsen_GEPCSAFTModel, V, T, z, n , _data=@f(data))
    _,_,_,_,_,m̄ = _data
    if n == 1
        corr = PCSAFTconsts.corr1
    elseif n == 2
        corr = PCSAFTconsts.corr2
    end
    res = zero(m̄)
    @inbounds for i ∈ 1:7
        corr1,corr2,corr3 = corr[i]
        ki = corr1 + (m̄-1)/m̄*corr2 + (m̄-1)/m̄*(m̄-2)/m̄*corr3
        res +=ki
    end
    return res
end

function Ii(model::Michelsen_GEPCSAFTModel, V, T, z, n , _data=@f(data))
    m = model.params.segment.values
    if n == 1
        corr = PCSAFTconsts.corr1
    elseif n == 2
        corr = PCSAFTconsts.corr2
    end
    res = zero(m)
    @inbounds for i ∈ 1:7
        corr1,corr2,corr3 = corr[i]
        ki = @. corr1 + (m-1)/m*corr2 + (m-1)/m*(m-2)/m*corr3
        res +=ki
    end
    return res
end

function newtons_method(f, df, xn)
    return xn - f/df 
end