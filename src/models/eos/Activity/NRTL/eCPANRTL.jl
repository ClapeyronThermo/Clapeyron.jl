struct eCPANRTLParam <: EoSParam
    dUref::PairParam{Float64}
    TU::PairParam{Float64}
    omegaU::PairParam{Float64}
    alpha::PairParam{Float64}
    b::SingleParam{Float64}
    Mw::SingleParam{Float64}
end

abstract type eCPANRTLModel <: ActivityModel end

struct eCPANRTL{c<:EoSModel} <: eCPANRTLModel
    components::Array{String,1}
    icomponents::UnitRange{Int}
    params::eCPANRTLParam
    puremodel::Vector{c}
    absolutetolerance::Float64
    references::Array{String,1}
end

@registermodel eCPANRTL

export eCPANRTL

function eCPANRTL(components::Vector{String}; puremodel=PR,
    userlocations=String[], 
     verbose=false)
    params = getparams(components, ["properties/critical.csv", "properties/molarmass.csv","Activity/eCPANRTL/eCPANRTL_unlike.csv","SAFT/CPA/CPA_like.csv"]; userlocations=userlocations, asymmetricparams=["dUref","TU","omegaU"], ignore_missing_singleparams=["dUref","TU","omegaU"], verbose=verbose)
    dUref  = params["dUref"]
    params["TU"].values[diagind(params["TU"].values)] .= 1.
    TU  = params["TU"]
    omegaU  = params["omegaU"]
    alpha = params["alpha"]
    Mw  = params["Mw"]
    b = params["b"]
    icomponents = 1:length(components)
    
    init_puremodel = [puremodel([components[i]]) for i in icomponents]
    packagedparams = eCPANRTLParam(dUref,TU,omegaU,alpha,b,Mw)
    references = String[]
    model = eCPANRTL(components,icomponents,packagedparams,init_puremodel,1e-12,references)
    return model
end

function excess_gibbs_free_energy(model::eCPANRTLModel,p,T,z)
    ΔUref = model.params.dUref.values
    TΔU = model.params.TU.values
    ωΔU = model.params.omegaU.values
    α = model.params.alpha.values
    b = model.params.b.values

    x = z ./ sum(z)

    τ = @. (ΔUref+ωΔU*((1-T/TΔU)^2-(1-298.15/TΔU)^2))/T
    G = @. exp(-α*τ)
    return sum(x[i]*sum(x[j]*b[j]*G[j,i]*τ[j,i] for j ∈ @comps)/sum(x[j]*b[j]*G[j,i] for j ∈ @comps) for i ∈ @comps)
end