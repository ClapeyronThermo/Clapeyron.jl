
abstract type AdvGEPCSAFTModel <: PCSAFTModel end


struct AdvGEPCSAFT{I <: IdealModel,T,γ} <: AdvGEPCSAFTModel
    components::Array{String,1}
    sites::SiteParam
    activity::γ
    params::PCSAFTParam{T}
    idealmodel::I
    assoc_options::AssocOptions
    references::Array{String,1}
end

"""
    AdvGEPCSAFT <: SAFTModel

    AdvGEPCSAFT(components;
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

Perturbed-Chain SAFT (PC-SAFT), with Gᴱ mixing rule - using the Michelsen (0 pressure) limit.

"""
AdvGEPCSAFT

export AdvGEPCSAFT
function AdvGEPCSAFT(components;
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

    packagedparams = PCSAFTParam(Mw, segment, sigma, epsilon, epsilon_assoc, bondvol)

    init_idealmodel = init_model(idealmodel,components,ideal_userlocations,verbose)
    init_activity = init_model(activity,components,activity_userlocations,verbose)
    references = String["10.1021/acs.iecr.2c03464"]
    model = AdvGEPCSAFT(format_components(components),sites,init_activity,packagedparams,init_idealmodel,assoc_options,references)
    set_reference_state!(model,reference_state;verbose)
    return model
end

function _pcsaft(model::AdvGEPCSAFT{I,T}) where {I,T}
    return PCSAFT{I,T}(model.components,model.sites,model.params,model.idealmodel,model.assoc_options,model.references)
end

function a_res(model::AdvGEPCSAFTModel, V, T, z)    
    _data = @f(data)
    return @f(a_hc,_data) + @f(a_disp,_data) + @f(a_assoc,_data)
end

function data(model::AdvGEPCSAFTModel,V,T,z)
    _d = @f(d)
    ζ0,ζ1,ζ2,ζ3 = @f(ζ0123,_d)
    m = model.params.segment.values
    m̄ = dot(z, m)/sum(z)
    return (_d,ζ0,ζ1,ζ2,ζ3,m̄)
end


# function a_hc(model::AdvGEPCSAFTModel, V, T, z , _data = @f(data))
#     _,_,_,_,η,m̄ = _data
#     g_hs = (1-η/2)/(1-η)^3
#     a_hs = (4η-3η^2)/(1-η)^2
#     return m̄*a_hs - (m̄-1)*log(g_hs)
# end

# function g_hs(model::AdvGEPCSAFTModel, V, T, z,_data = @f(data))
#     _,_,_,_,η,_ = _data
#     return (1-η/2)/(1-η)^3
# end

# function a_hs(model::AdvGEPCSAFTModel, V, T, z)
#     _,_,_,_,η,_ = _data
#     return (4η-3η^2)/(1-η)^2
# end


function m2ϵσ3(model::AdvGEPCSAFTModel, V, T, z, _data=@f(data))

    function q_i(α, m)
        c = [-1, -3.312911261372218, (log(m) + -2.26053639783524) * (4.373514008544666 - log(m)^2 * -0.5438472133955807 - 2.7803145680981665*log(m))]

        # c = [2.4943621118539628*(log(b))^2 + 317.1749262783832*log(b) + 10067.759452498541, 8.066923060464152*(log(b))^2 + 1065.837604157669*log(b) + 35238.98020488654]
        # return c[1]*α + c[2]
        return c[1]*α^2 + c[2]*α + c[3]
    end

    function α_mix(q̄,m̄)
        c = [-1, -3.312911261372218, (log(m̄) + -2.26053639783524) * (4.373514008544666 - log(m̄)^2 * -0.5438472133955807 - 2.7803145680981665*log(m̄))- q̄]
        # c = [2.4943621118539628*(log(b̄))^2 + 317.1749262783832*log(b̄) + 10067.759452498541, 8.066923060464152*(log(b̄))^2 + 1065.837604157669*log(b̄) + 35238.98020488654]
        
        # We have to solve the equation q(α, b) = q̄ 
        # where q(α, b) = c[1]*α + c[2]
        # return (q̄ - c[2])/c[1]
        return (-c[2]-sqrt(c[2]^2-4*c[1]*c[3]))/(2c[1])
    end

    di,ζ0,ζ1,ζ2,ζ3,m̄ = _data
    Tnum = promote_type(eltype(z), typeof(V), typeof(T))
    act_model = typeof(model.activity).name.wrapper
    N = length(z)

    m = model.params.segment.values
    ϵ = diagvalues(model.params.epsilon)
    σ = diagvalues(model.params.sigma)
    α = ϵ./T

    b = m.*di.^3

    q = @. q_i(α, m)*m
    # println(q)
    # println(α)
    # println(b)
    # println(m)

    Σz = sum(z)
    #Iᵢ = @f(Ii,1,_data)
    b̄ = zero(Base.promote_eltype(model,V,T,z))
    m²σ³ = zero(b̄)
    A = zero(b̄)
    B = zero(b̄)
    @inbounds for i ∈ @comps
        mᵢ,bᵢ,σᵢ,qᵢ,zᵢ = m[i],b[i],σ[i],q[i],z[i]
        σ³ᵢ = σᵢ*σᵢ*σᵢ
        m²σ³ += zᵢ*mᵢ*σ³ᵢ
        b̄ += zᵢ*bᵢ
        A += zᵢ*qᵢ
        B += zᵢ*log(bᵢ)
    end
    m²σ³,b̄ = m²σ³/Σz,b̄/Σz
    A, B = A/Σz, B/Σz
    gₑ = excess_gibbs_free_energy(model.activity,V,T,z)/(R̄*T*Σz)
    
    q̄ = gₑ + (log(b̄) - B) +  A
    ᾱ = α_mix(q̄/m̄, m̄)
    # println("g_E/RT = ", gₑ)
    # println("log( b̄ ) = ", log(b̄))
    # println("B = ", B)
    # println("A = ", A)
    # println("q̄ = ", q̄)
    # println("ᾱ = ", ᾱ)
    m2ϵσ3₁ = ᾱ*m²σ³*m̄
    m2ϵσ3₂ = ᾱ*m2ϵσ3₁

    return m2ϵσ3₁, m2ϵσ3₂
end