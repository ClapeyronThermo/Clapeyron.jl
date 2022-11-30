struct eCPANRTLParam <: EoSParam
    A::PairParam{Float64}
    B::PairParam{Float64}
    C::PairParam{Float64}
    alpha::PairParam{Float64}
    a::PairParam{Float64}
    b::PairParam{Float64}
    c1::SingleParam{Float64}
    Tc::SingleParam{Float64}
    Mw::SingleParam{Float64}
end

abstract type eCPANRTLModel <: ActivityModel end

struct eCPANRTL{c<:EoSModel} <: eCPANRTLModel
    components::Array{String,1}
    params::eCPANRTLParam
    puremodel::Vector{c}
    absolutetolerance::Float64
    references::Array{String,1}
end

@registermodel eCPANRTL

export eCPANRTL

function eCPANRTL(components::Vector{String}; puremodel=PR,
    userlocations=String[],
    pure_userlocations = String[], 
    verbose=false)

    params = getparams(components, ["properties/critical.csv", "properties/molarmass.csv","Activity/eCPANRTL/eCPANRTL_unlike.csv","SAFT/sCPA/sCPA_like.csv","SAFT/sCPA/sCPA_unlike.csv"]; userlocations=userlocations, asymmetricparams=["A","B","C"], ignore_missing_singleparams=["A","B","C"], verbose=verbose)
    A  = params["A"]
    B  = params["B"]
    C  = params["C"]
    alpha = params["alpha"]
    Mw  = params["Mw"]
    k  = params["k"]
    Tc = params["Tc"]
    c1 = params["c1"]
    params["a"].values .*= 1E-1
    params["b"].values .*= 1E-3
    a  = epsilon_LorentzBerthelot(params["a"], k)
    b  = sigma_LorentzBerthelot(params["b"])
    
    
    pure = init_puremodel(puremodel,components,pure_userlocations,verbose)
    packagedparams = eCPANRTLParam(A,B,C,alpha,a,b,c1,Tc,Mw)
    references = String[]
    model = eCPANRTL(components,packagedparams,pure,1e-12,references)
    return model
end

function excess_gibbs_free_energy(model::eCPANRTLModel,p,T,z)
    A = model.params.A.values
    B = model.params.B.values
    C = model.params.C.values
    α = model.params.alpha.values
    b = model.params.b.diagvalues
    a₀ = model.params.a.values
    c1 = model.params.c1.values
    Tc = model.params.Tc.values

    Tr = @. T/Tc

    ᾱ = @. (1+c1*(1-√(Tr)))

    a = sqrt.(a₀.*a₀').*(ᾱ.*ᾱ')*log(2)/R̄

    ΔU = sign.(diag(a)).*((diag(a)./b)'.-a./b)

    x = z ./ sum(z)

    τ = @. (A+B*T+C*T^2+ΔU)/T
    G = @. exp(-α*τ)
    return sum(x[i]*sum(x[j]*b[j]*G[j,i]*τ[j,i] for j ∈ @comps)/sum(x[j]*b[j]*G[j,i] for j ∈ @comps) for i ∈ @comps)*sum(z)*R̄*T
end