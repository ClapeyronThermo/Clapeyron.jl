abstract type MHV1RuleModel <: ActivityMixingRule end

struct MHV1Rule{γ} <: MHV1RuleModel
    components::Array{String,1}
    activity::γ
    references::Array{String,1}
end

"""
    MHV1Rule{γ} <: MHV1RuleModel

    MHV1Rule(components;
    activity = Wilson,
    userlocations = String[],
    activity_userlocations = String[],
    verbose::Bool=false)

## Input Parameters

None

## Input models

- `activity`: Activity Model

## Description

Modified Huron-Vidal Mixing Rule, First Order
```
bᵢⱼ = (1 - lᵢⱼ)(bᵢ + bⱼ)/2
b̄ = ∑bᵢⱼxᵢxⱼ
c̄ = ∑cᵢxᵢ
ā = b̄RT(∑[xᵢaᵢᵢαᵢ/(RTbᵢᵢ)] - [gᴱ/RT + ∑log(bᵢᵢ/b̄)]/q)

if the model is Peng-Robinson:
    q = 0.53
if the model is Redlich-Kwong:
    q = 0.593
```

to use different values for `q`, overload `Clapeyron.MHV1q(::CubicModel,::MHV1Model) = q`

## Model Construction Examples
```
# Using the default database
mixing = MHV1Rule(["water","carbon dioxide"]) #default: Wilson Activity Coefficient.
mixing = MHV1Rule(["water","carbon dioxide"],activity = NRTL) #passing another Activity Coefficient Model.
mixing = MHV1Rule([("ethane",["CH3" => 2]),("butane",["CH2" => 2,"CH3" => 2])],activity = UNIFAC) #passing a GC Activity Coefficient Model.

# Passing a prebuilt model

act_model = NRTL(["water","ethanol"],userlocations = (a = [0.0 3.458; -0.801 0.0],b = [0.0 -586.1; 246.2 0.0], c = [0.0 0.3; 0.3 0.0]))
mixing = MHV1Rule(["water","ethanol"],activity = act_model)

# Using user-provided parameters

# Passing files or folders
mixing = MHV1Rule(["water","ethanol"]; activity = NRTL, activity_userlocations = ["path/to/my/db","nrtl_ge.csv"])

# Passing parameters directly
mixing = MHV1Rule(["water","ethanol"];
                activity = NRTL,
                userlocations = (a = [0.0 3.458; -0.801 0.0],
                    b = [0.0 -586.1; 246.2 0.0],
                    c = [0.0 0.3; 0.3 0.0])
                )
```

## References
1. Michelsen, M. L. (1990). A modified Huron-Vidal mixing rule for cubic equations of state. Fluid Phase Equilibria, 60(1–2), 213–219. [doi:10.1016/0378-3812(90)85053-d](https://doi.org/10.1016/0378-3812(90)85053-d)

"""
MHV1Rule

export MHV1Rule
function MHV1Rule(components; activity = Wilson, userlocations = String[],activity_userlocations = String[], verbose::Bool=false)
    _activity = init_mixing_act(activity,components,activity_userlocations,verbose)
    references = ["10.1016/0378-3812(90)85053-D"]
    model = MHV1Rule(format_components(components), _activity,references)
    return model
end

MHV1q(::MHV1RuleModel,::PRModel) = 0.53
MHV1q(::MHV1RuleModel,::RKModel) = 0.593 #check if it applies to vanilla RK
#MHV1q(::MHV1RuleModel,::SRKModel) = 0.593 causes ambiguities
MHV1q(::MHV1RuleModel,::vdWModel) = 0.85

function mixing_rule(model::ABCubicModel,V,T,z,mixing_model::MHV1RuleModel,α,a,b)
    n = sum(z)
    invn = (one(n)/n)
    invn2 = invn^2
    g_E = excess_gibbs_free_energy(mixing_model.activity,1e5,T,z) / n
    b̄ = dot(z,Symmetric(b),z) * invn2
    c̄ = translation2(model,V,T,z,model.translation,a,b,α)*invn
    q = MHV1q(mixing_model,model)
    Σlogb = sum(z[i]*log(b̄/b[i,i]) for i ∈ @comps)*invn
    Σab = invn*sum(z[i]*a[i,i]*α[i]/b[i,i]/(R̄*T) for i ∈ @comps)
    ā = b̄*R̄*T*(Σab-1/q*(g_E/(R̄*T)+Σlogb))
    return ā,b̄,c̄
end

function cubic_get_l(model::CubicModel,mixing::MHV1RuleModel,params)
    return get_k_mean(params.b.values)
end
