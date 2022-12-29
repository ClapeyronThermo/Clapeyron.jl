abstract type MHV2RuleModel <: ActivityMixingRule end

struct MHV2Rule{γ} <: MHV2RuleModel
    components::Array{String,1}
    activity::γ
    references::Array{String,1}
end

@registermodel MHV2Rule

"""
    MHV2Rule{γ} <: MHV2RuleModel

    MHV2Rule(components::Vector{String};
    activity = Wilson,
    userlocations::Vector{String}=String[],
    activity_userlocations::Vector{String}=String[],
    verbose::Bool=false)

## Input Parameters

None

## Input models

- `activity`: Activity Model

## Description

Modified Huron-Vidal Mixing Rule, Second Order.
```
aᵢⱼ = √(aᵢaⱼ)(1 - kᵢⱼ)
bᵢⱼ = (bᵢ + bⱼ)/2
b̄ = ∑bᵢⱼxᵢxⱼ
c̄ = ∑cᵢxᵢ
ᾱᵢ  = aᵢαᵢ/bᵢRT
ċ = -q₁*Σᾱᵢxᵢ - q₂*Σᾱᵢxᵢ^2 - gᴱ/RT - ∑log(bᵢᵢ/b̄)
ā = (-q₁ - √(q₁^2 - 4q₂ċ))/(2q₂)

if the model is Peng-Robinson:
    q₁ = -0.4347, q₂ = -0.003654
if the model is Redlich-Kwong:
    q₁ = -0.4783, q₂ = -0.0047
    (-0.4783,-0.0047)
```

to use different values for `q₁` and `q₂`, overload `Clapeyron.MHV1q(::CubicModel,::MHV2Model) = (q₁,q₂)`


## References
1. Michelsen, M. L. (1990). A modified Huron-Vidal mixing rule for cubic equations of state. Fluid Phase Equilibria, 60(1–2), 213–219. [doi:10.1016/0378-3812(90)85053-d](https://doi.org/10.1016/0378-3812(90)85053-d)

"""
MHV2Rule


export MHV2Rule
function MHV2Rule(components::Vector{String}; activity = Wilson, userlocations::Vector{String}=String[],activity_userlocations::Vector{String}=String[], verbose::Bool=false)
    _activity = init_model(activity,components,activity_userlocations,verbose)

    references = ["10.1016/0378-3812(90)85053-D"]
    model = MHV2Rule(components, _activity,references)
    return model
end

MHV2q(::MHV2RuleModel,::PRModel) = (-0.4347,-0.003654)
MHV2q(::MHV2RuleModel,::RKModel) = (-0.4783,-0.0047)

function mixing_rule(model::Union{PRModel,RKModel},V,T,z,mixing_model::MHV2RuleModel,α,a,b,c)
    n = sum(z)
    invn = (one(n)/n)
    invn2 = invn^2
    g_E = excess_gibbs_free_energy(mixing_model.activity,1e5,T,z)*invn
    b̄ = dot(z,Symmetric(b),z) * invn2
    c̄ = dot(z,c)*invn
    #ᾱ = a.*sqrt.(α.*α')./(b*R̄*T)
    q1,q2 = MHV2q(mixing_model,model)
    Σlogb = zero(first(z))
    Σab = zero(T+first(z))
    Σab2 = Σab
    for i ∈ @comps
        ᾱ  = a[i,i]*α[i]/(b[i,i]*R̄*T)
        zia = z[i]*ᾱ
        Σab += zia
        Σab2 += zia*ᾱ
        Σlogb += z[i]*log(b̄/b[i,i])
    end
    Σab *= invn
    Σab2 *= invn
    Σlogb *= invn
    c  = -q1*Σab-q2*Σab2-g_E/(R̄*T)-Σlogb
    ā = b̄*R̄*T*(-q1-sqrt(q1^2-4*q2*c))/(2*q2)
    return ā,b̄,c̄
end
#=
function mixing_rule(model::PRModel,V,T,z,mixing_model::MHV2RuleModel,α,a,b,c)
    n = sum(z)
    x = z./n
    invn2 = (one(n)/n)^2
    g_E = excess_gibbs_free_energy(mixing_model.activity,1e5,T,z) / n
    b̄ = dot(z,Symmetric(b),z) * invn2
    c̄ = dot(z,c)/n

    ᾱ = a.*sqrt.(α.*α')./(b*R̄*T)

    q1 = -0.4347
    q2 = -0.003654
    c  = -q1*sum(x[i]*ᾱ[i,i] for i ∈ @comps)-q2*sum(x[i]*ᾱ[i,i]^2 for i ∈ @comps)-g_E/(R̄*T)-sum(x[i]*log(b̄/b[i,i]) for i ∈ @comps)

    ā = b̄*R̄*T*(-q1-sqrt(q1^2-4*q2*c))/(2*q2)
    return ā,b̄,c̄
end=#
