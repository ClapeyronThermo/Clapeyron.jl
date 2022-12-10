abstract type MHV1RuleModel <: ActivityMixingRule end

struct MHV1Rule{γ} <: MHV1RuleModel
    components::Array{String,1}
    activity::γ
    references::Array{String,1}
end

@registermodel MHV1Rule

"""
    MHV1Rule{γ} <: MHV1RuleModel

    MHV1Rule(components::Vector{String};
    activity = Wilson,
    userlocations::Vector{String}=String[],
    activity_userlocations::Vector{String}=String[],
    verbose::Bool=false)

## Input Parameters

None

## Input models

- `activity`: Activity Model

## Description

Modified Huron-Vidal Mixing Rule, First Order
```
aᵢⱼ = √(aᵢaⱼ)(1-kᵢⱼ)
bᵢⱼ = (bᵢ + bⱼ)/2
b̄ = ∑bᵢⱼxᵢxⱼ
c̄ = ∑cᵢxᵢ
ā = b̄RT(∑[xᵢaᵢᵢαᵢ/(RTbᵢᵢ)] - [gᴱ/RT + ∑log(bᵢᵢ/b̄)]/q)

if the model is Peng-Robinson:
    q = 0.53
if the model is Redlich-Kwong:
    q = 0.593
```

to use different values for `q`, overload `Clapeyron.MHV1q(::CubicModel,::MHV1Model) = q`


## References
1. Michelsen, M. L. (1990). A modified Huron-Vidal mixing rule for cubic equations of state. Fluid Phase Equilibria, 60(1–2), 213–219. [doi:10.1016/0378-3812(90)85053-d](https://doi.org/10.1016/0378-3812(90)85053-d)

"""
MHV1Rule

export MHV1Rule
function MHV1Rule(components::Vector{String}; activity = Wilson, userlocations::Vector{String}=String[],activity_userlocations::Vector{String}=String[], verbose::Bool=false)
    _activity = init_model(activity,components,activity_userlocations,verbose)
    references = ["10.1016/0378-3812(90)85053-D"]
    model = MHV1Rule(components, _activity,references)
    return model
end

MHV1q(::MHV1RuleModel,::PRModel) = 0.53
MHV1q(::MHV1RuleModel,::RKModel) = 0.593 #check if it applies to vanilla RK
#MHV1q(::MHV1RuleModel,::SRKModel) = 0.593 causes ambiguities
MHV1q(::MHV1RuleModel,::vdWModel) = 0.85

function mixing_rule(model::ABCubicModel,V,T,z,mixing_model::MHV1RuleModel,α,a,b,c)
    n = sum(z)
    invn = (one(n)/n)
    invn2 = invn^2
    g_E = excess_gibbs_free_energy(mixing_model.activity,1e5,T,z) / n
    b̄ = dot(z,Symmetric(b),z) * invn2
    c̄ = dot(z,c)/n
    q = MHV1q(mixing_model,model)
    Σlogb = sum(z[i]*log(b̄/b[i,i]) for i ∈ @comps)*invn
    Σab = invn*sum(z[i]*a[i,i]*α[i]/b[i,i]/(R̄*T) for i ∈ @comps)
    ā = b̄*R̄*T*(Σab-1/q*(g_E/(R̄*T)+Σlogb))
    return ā,b̄,c̄
end
