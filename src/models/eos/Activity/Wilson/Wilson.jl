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

has_sites(::Type{<:WilsonModel}) = false
has_groups(::Type{<:WilsonModel}) = false
built_by_macro(::Type{<:WilsonModel}) = false

function Base.show(io::IO, mime::MIME"text/plain", model::Wilson)
    return eosshow(io, mime, model)
end

function Base.show(io::IO, model::Wilson)
    return eosshow(io, model)
end

Base.length(model::Wilson) = Base.length(model.icomponents)

molecular_weight(model::Wilson,z=SA[1.0]) = comp_molecular_weight(mw(model),z)

export Wilson

function Wilson(components::Vector{String}; puremodel=PR,
    userlocations=String[], 
     verbose=false)
    params = getparams(components, ["properties/critical.csv", "properties/molarmass.csv","Activity/Wilson/Wilson_unlike.csv"]; userlocations=userlocations, asymmetricparams=["g"], ignore_missing_singleparams=true, verbose=verbose)
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
    V =  @. (R̄ *Tc/Pc)*pow(ZRA,(1 + pow((1-Tr),2/7)))
    Λ = (V' ./ V) .*exp.(-model.params.g.values/R̄/T)
    x = z ./ sum(z)
    lnγ = 1 .- log.(sum(x[i]*Λ[:,i] for i ∈ @comps)) .-sum(x[j] .*Λ[j,:] ./(sum(x[i]*Λ[j,i] for i ∈ @comps)) for j ∈ @comps)
    return exp.(lnγ)
end

function excess_gibbs_free_energy(model::WilsonModel,p,T,z)
    ZRA = model.params.ZRA.values
    Tc  = model.params.Tc.values
    Pc  = model.params.Pc.values
    
    Tr  = T ./ Tc
    V =  @. (R̄ *Tc/Pc)*pow(ZRA,(1 + pow((1-Tr),2/7)))
    Λ = @. (V' / V) *exp(-model.params.g.values/R̄/T)
    x = z ./ sum(z)
    g_E = -sum(x[j]*log(sum(x[i]*Λ[j,i] for i ∈ @comps)) for j ∈ @comps)
    return g_E*sum(z)*R̄*T
end

function bubble_pressure(model::ActivityModel,T,x)
    sat = sat_pure.(model.puremodel,T)
    p_sat = [tup[1] for tup in sat]
    γ     = activity_coefficient(model,1e-4,T,x)
    p     = sum(x.*γ.*p_sat)
    y     = x.*γ.*p_sat ./ p
    return (p,y)
end

is_splittable(::Wilson) = true