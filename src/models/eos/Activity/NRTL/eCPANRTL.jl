struct eCPANRTLParam <: EoSParam
    A::PairParam{Float64}
    B::PairParam{Float64}
    C::PairParam{Float64}
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
    params = getparams(components, ["properties/critical.csv", "properties/molarmass.csv","Activity/eCPANRTL/eCPANRTL_unlike.csv","SAFT/sCPA/sCPA_like.csv"]; userlocations=userlocations, asymmetricparams=["A","B","C"], ignore_missing_singleparams=["A","B","C"], verbose=verbose)
    A  = params["A"]
    B  = params["B"]
    C  = params["C"]
    alpha = params["alpha"]
    Mw  = params["Mw"]
    b = params["b"]
    icomponents = 1:length(components)
    
    init_puremodel = [puremodel([components[i]]) for i in icomponents]
    packagedparams = eCPANRTLParam(A,B,C,alpha,b,Mw)
    references = String[]
    model = eCPANRTL(components,icomponents,packagedparams,init_puremodel,1e-12,references)
    return model
end

function excess_gibbs_free_energy(model::eCPANRTLModel,p,T,z)
    A = model.params.A.values
    B = model.params.B.values
    C = model.params.C.values
    α = model.params.alpha.values
    b = model.params.b.values

    x = z ./ sum(z)

    τ = @. (A+B*T+C*T^2)/T
    G = @. exp(-α*τ)
    return sum(x[i]*sum(x[j]*b[j]*G[j,i]*τ[j,i] for j ∈ @comps)/sum(x[j]*b[j]*G[j,i] for j ∈ @comps) for i ∈ @comps)*sum(z)*R̄*T
end