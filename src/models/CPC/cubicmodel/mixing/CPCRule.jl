
struct CPCRuleParam  <: EoSParam
    segment::SingleParam{Float64}
    omega_a::SingleParam{Float64}
    omega_b::SingleParam{Float64}
end

struct CPCRule{R} <: CPCRuleModel
    components::Vector{String}
    params::CPCRuleParam
    rdf::R
    references::Vector{String}
end

@registermodel CPCRule

"""
    CPCRule <: CPCRuleModel
    
    CPCRule(components::Vector{String};
    userlocations::Vector{String}=String[],
    rdf = ElliottRDF()
    verbose::Bool=false)

## Input Parameters

None

## Input models 

- `rdf`: CPC' Radial distance function model

## Description

Mixing Rule used by Cubic plus Chain equation of state.
```
aᵢⱼ = √(aᵢaⱼ)(1-kᵢⱼ)
bᵢⱼ = (bᵢ + bⱼ)/2
m̄b̄ = ∑bᵢxᵢmᵢ
m̄²ā = ∑aᵢⱼxᵢxⱼmᵢmⱼ
c̄ = ∑cᵢxᵢ
```
"""
CPCRule
export CPCRule
function CPCRule(components::Vector{String}; userlocations::Vector{String}=String[],rdf = ElliottRDF(), verbose::Bool=false)
    params = getparams(components,["todo"])
    segment = params["m"]
    references = ["10.1021/acs.iecr.9b00435"]
    omega_a = SingleParam("Ωa",components)
    omega_b = SingleParam("Ωb",components)
    params = CPCRuleParam(segment,omega_a,omega_b)
    model = CPCRule(components, params, rdf, references)
    return model
end



function mixing_rule(model::ABCubicModel,V,T,z,mixing_model::CPCRuleModel,α,a,b,c)
    n = sum(z)
    invn = (one(n)/n)
    invn2 = invn^2
    m = mixing_model.params.segment
    invn = (one(n)/n)
    invn2 = invn^2
    #b̄ = dot(z,Symmetric(b),z) * invn2
    ā,b̄ = zero(T+first(z)),zero(first(z))
    m̄,c̄ = dot(m,z)*invn,dot(z,c)*invn
    for i in 1:length(z)
        zi,αi,mi = z[i],α[i],m[i]
        zi2 = zi^2
        b̄ += mi*b[i,i]*zi2
        ā += mi*mia[i,i]*αi*zi^2
        for j in 1:(i-1)
            zij = zi*z[j]
            ā += 2*mi*m[j]*a[i,j]*sqrt(αi*α[j])*zij
        end
    end
    ā *= invn2/(m̄*m̄)
    b̄ *= invn2/m̄
    #dot(z,Symmetric(a .* sqrt.(α*α')),z) * invn2
    return ā,b̄,c̄
end