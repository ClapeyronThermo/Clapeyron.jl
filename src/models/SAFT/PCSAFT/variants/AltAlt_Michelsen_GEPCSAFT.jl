
abstract type AltAltAdvGEPCSAFTModel <: PCSAFTModel end


struct AltAltAdvGEPCSAFT{I <: IdealModel,T,γ} <: AltAltAdvGEPCSAFTModel
    components::Array{String,1}
    sites::SiteParam
    activity::γ
    params::PCSAFTParam{T}
    idealmodel::I
    assoc_options::AssocOptions
    references::Array{String,1}
    λ::T
end

"""
    AltAltAdvGEPCSAFT <: SAFTModel

    AltAltAdvGEPCSAFT(components;
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
AltAltAdvGEPCSAFT

export AltAltAdvGEPCSAFT
function AltAltAdvGEPCSAFT(components;
    idealmodel = BasicIdeal,
    activity = UNIFAC,
    userlocations = String[],
    ideal_userlocations = String[],
    activity_userlocations = String[],
    assoc_options = AssocOptions(),
    reference_state = nothing,
    λ = 1.0,
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
    model = AltAltAdvGEPCSAFT(format_components(components),sites,init_activity,packagedparams,init_idealmodel,assoc_options,references,λ)
    set_reference_state!(model,reference_state;verbose)
    return model
end

function _pcsaft(model::AltAltAdvGEPCSAFT{I,T}) where {I,T}
    return PCSAFT{I,T}(model.components,model.sites,model.params,model.idealmodel,model.assoc_options,model.references)
end

function m2ϵσ3(model::AltAltAdvGEPCSAFTModel, V, T, z, _data=@f(data))

    function q_i(α, b)
        c = [1, 0.15498301934788844, 13.13195034750734, 1.575271648063568, 2, -12.509610376473578]
        (c[1] + c[2]*log(b))*α^2 + (c[3] + c[4]*log(b))*α + c[5]*log(b) + c[6]
    end

    function α_mix(q̄,b̄)
        c = [1, 0.15498301934788844, 13.13195034750734, 1.575271648063568, 2, -12.509610376473578]
        A = (c[1] + c[2]*log(b̄))
        B = (c[3] + c[4]*log(b̄))
        C = c[5]*log(b̄) + c[6] - q̄
        # Solve the quadratic equation A*α^2 + B*α + C = 0
        return (-B-sqrt(B^2 - 4*A*C))/(2*A)
    end

    di,ζ0,ζ1,ζ2,ζ3,m̄ = _data
    Tnum = promote_type(eltype(z), typeof(V), typeof(T))
    act_model = typeof(model.activity).name.wrapper
    N = length(z)

    m = model.params.segment.values
    ϵ = diagvalues(model.params.epsilon)
    σ = diagvalues(model.params.sigma)
    α = m.*ϵ./T

    b = m.*di.^3

    q = @. q_i(α, b)
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
        m²σ³ += zᵢ*mᵢ*mᵢ*σ³ᵢ
        b̄ += zᵢ*bᵢ
        A += zᵢ*qᵢ
        B += zᵢ*log(bᵢ)
    end
    m²σ³,b̄ = m²σ³/Σz,b̄/Σz
    A, B = A/Σz, B/Σz
    gₑ = excess_gibbs_free_energy(model.activity,V,T,z)/(R̄*T*Σz)
    q̄ = gₑ + model.λ*(log(b̄)-B) +  A
    ᾱ = α_mix(q̄, b̄)
    m2ϵσ3₁ = ᾱ*m²σ³/m̄
    m2ϵσ3₂ = m2ϵσ3₁*m2ϵσ3₁/m²σ³

    return m2ϵσ3₁, m2ϵσ3₂
end