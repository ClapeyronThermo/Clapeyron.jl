struct WilsonParam <: EoSParam
    g::PairParam{Float64}
    Tc::SingleParam{Float64}
    Pc::SingleParam{Float64}
    ZRA::SingleParam{Float64}
    Mw::SingleParam{Float64}
end

abstract type WilsonModel <: ActivityModel end

struct Wilson{c<:EoSModel} <: WilsonModel
    components::Array{String,1}
    icomponents::UnitRange{Int}
    params::WilsonParam
    puremodel::Vector{c}
    absolutetolerance::Float64
    references::Array{String,1}
end

@registermodel Wilson
export Wilson

function Wilson(components::Vector{String}; puremodel=PR,
    userlocations=String[], 
     verbose=false)
    params = getparams(components, ["properties/critical.csv", "properties/molarmass.csv","Activity/Wilson/Wilson_unlike.csv"]; userlocations=userlocations, asymmetricparams=["g"], ignore_missing_singleparams=["g"], verbose=verbose)
    g  = params["g"]
    Tc        = params["Tc"]
    pc        = params["pc"]
    Mw        = params["Mw"]
    ZRA       = SingleParam(params["w"],"acentric factor")
    ZRA.values .*= -0.08775
    ZRA.values .+= 0.29056
    icomponents = 1:length(components)
    
    init_puremodel = [puremodel([components[i]]) for i in icomponents]
    packagedparams = WilsonParam(g,Tc,pc,ZRA,Mw)
    references = String[]
    model = Wilson(components,icomponents,packagedparams,init_puremodel,1e-12,references)
    return model
end

function activity_coefficient(model::WilsonModel,p,T,z)
    ZRA = model.params.ZRA.values
    Tc  = model.params.Tc.values
    Pc  = model.params.Pc.values
    
    Tr  = T ./ Tc
    V =  @. (R̄ *Tc/Pc)*ZRA^(1 + (1-Tr)^2/7)
    Λ = (V' ./ V) .*exp.(-model.params.g.values/R̄/T)
    x = z ./ sum(z)
    lnγ = 1 .- log.(sum(x[i]*Λ[:,i] for i ∈ @comps)) .-sum(x[j] .*Λ[j,:] ./(sum(x[i]*Λ[j,i] for i ∈ @comps)) for j ∈ @comps)
    return exp.(lnγ)
end