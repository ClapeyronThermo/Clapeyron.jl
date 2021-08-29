struct NRTLParam <: EoSParam
    a::PairParam{Float64}
    b::PairParam{Float64}
    c::PairParam{Float64}
    Mw::SingleParam{Float64}
end

abstract type NRTLModel <: ActivityModel end

struct NRTL{c<:EoSModel} <: NRTLModel
    components::Array{String,1}
    icomponents::UnitRange{Int}
    params::NRTLParam
    puremodel::Vector{c}
    absolutetolerance::Float64
    references::Array{String,1}
end

has_sites(::Type{<:NRTLModel}) = false
has_groups(::Type{<:NRTLModel}) = false
built_by_macro(::Type{<:NRTLModel}) = false

function Base.show(io::IO, mime::MIME"text/plain", model::NRTL)
    return eosshow(io, mime, model)
end

function Base.show(io::IO, model::NRTL)
    return eosshow(io, model)
end

Base.length(model::NRTL) = Base.length(model.icomponents)

molecular_weight(model::NRTL,z=SA[1.0]) = comp_molecular_weight(mw(model),z)

export NRTL

function NRTL(components::Vector{String}; puremodel=PR,
    userlocations=String[], 
     verbose=false)
    params = getparams(components, ["properties/critical.csv", "properties/molarmass.csv","Activity/NRTL/NRTL_unlike.csv"]; userlocations=userlocations, asymmetricparams=["a","b"], ignore_missing_singleparams=["a","b"], verbose=verbose)
    a  = params["a"]
    b  = params["b"]
    c  = params["c"]
    Mw  = params["Mw"]
    icomponents = 1:length(components)
    
    init_puremodel = [puremodel([components[i]]) for i in icomponents]
    packagedparams = NRTLParam(a,b,c,Mw)
    references = String[]
    model = NRTL(components,icomponents,packagedparams,init_puremodel,1e-12,references)
    return model
end

function activity_coefficient(model::NRTLModel,p,T,z)
    a = model.params.a.values
    b = model.params.b.values
    c = model.params.c.values

    x = z ./ sum(z)

    τ = @. a+b/T
    G = @. exp(-c*τ)
    lnγ = sum(x[j]*τ[j,:].*G[j,:] for j ∈ @comps)./sum(x[k]*G[k,:] for k ∈ @comps)+sum(x[j]*G[:,j]/sum(x[k]*G[k,j] for k ∈ @comps).*(τ[:,j] .-sum(x[m]*τ[m,j]*G[m,j] for m ∈ @comps)/sum(x[k]*G[k,j] for k ∈ @comps)) for j in @comps)
    return exp.(lnγ)
end

is_splittable(::NRTL) = true
