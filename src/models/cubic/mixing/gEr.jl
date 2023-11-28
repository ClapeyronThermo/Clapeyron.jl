abstract type gErRuleModel <: ActivityMixingRule end

struct gErRule{γ} <: gErRuleModel
    components::Array{String,1}
    activity::γ
    references::Array{String,1}
end

"""
    gErRule{γ} <: gErRuleModel

    gErRule(components;
    activity = NRTL,
    userlocations=String[],
    activity_userlocations=String[],
    verbose::Bool=false)

## Input Parameters

None

## Input models

- `activity`: Activity Model

## Description

Mixing rule that uses the residual part of the activity coefficient model:
```
bᵢⱼ = ( (bᵢ^(2/3) + bⱼ^(2/3)) / 2 )^(3/2)
b̄ = ∑bᵢⱼxᵢxⱼ
c̄ = ∑cᵢxᵢ
ā/b̄ = b̄RT*(∑xᵢaᵢ/bᵢ + gᴱᵣ/Λ
Λ = 1/(r₂ - r₁) * log((1 - r₂)/(1 - r₁))
```

## Model Construction Examples
```
# Using the default database
mixing = gErRule(["water","carbon dioxide"]) #default: NRTL Activity Coefficient.
mixing = gErRule(["water","carbon dioxide"],activity = Wilson) #passing another Activity Coefficient Model.
mixing = gErRule([("ethane",["CH3" => 2]),("butane",["CH2" => 2,"CH3" => 2])],activity = UNIFAC) #passing a GC Activity Coefficient Model.

# Passing a prebuilt model

act_model = NRTL(["water","ethanol"],userlocations = (a = [0.0 3.458; -0.801 0.0],b = [0.0 -586.1; 246.2 0.0], c = [0.0 0.3; 0.3 0.0]))
mixing = gErRule(["water","ethanol"],activity = act_model)

# Using user-provided parameters

# Passing files or folders
mixing = gErRule(["water","ethanol"]; activity = NRTL, activity_userlocations = ["path/to/my/db","nrtl_ge.csv"])

# Passing parameters directly
mixing = gErRule(["water","ethanol"];
                activity = NRTL,
                userlocations = (a = [0.0 3.458; -0.801 0.0],
                    b = [0.0 -586.1; 246.2 0.0],
                    c = [0.0 0.3; 0.3 0.0])
                )
```


## References
1. Piña-Martinez, A., Privat, R., Nikolaidis, I. K., Economou, I. G., & Jaubert, J.-N. (2021). What is the optimal activity coefficient model to be combined with the translated–consistent Peng–Robinson equation of state through advanced mixing rules? Cross-comparison and grading of the Wilson, UNIQUAC, and NRTL aE models against a benchmark database involving 200 binary systems. Industrial & Engineering Chemistry Research, 60(47), 17228–17247. [doi:10.1021/acs.iecr.1c03003](https://doi.org/10.1021/acs.iecr.1c03003)
"""
gErRule

default_references(::Type{gErRule}) =  ["10.1016/0378-3812(90)85053-D"]

export gErRule
function gErRule(components; activity = NRTL, userlocations=String[],activity_userlocations=String[], verbose::Bool=false)
    _components = format_components(components)
    _activity = init_model(activity,_components,activity_userlocations,verbose)
    model = gErRule(_components, _activity,default_references(gErRule))
    return model
end

function ab_premixing(model::PRModel,mixing::gErRuleModel,k = nothing, l = nothing)
    Ωa, Ωb = ab_consts(model)
    _Tc = model.params.Tc
    _pc = model.params.Pc
    a = model.params.a
    b = model.params.b
    diagvalues(a) .= @. Ωa*R̄^2*_Tc^2/_pc
    diagvalues(b) .= @. Ωb*R̄*_Tc/_pc
    epsilon_LorentzBerthelot!(a)
    gEr_mix(bi,bj,lij) = mix_powmean(bi,bj,lij,2/3)
    kij_mix!(gEr_mix,b,l)
    return a,b
end

function cubic_get_l(model::CubicModel,mixing::gErRuleModel,params)
    return get_k_powmean(params.b.values,2/3)
end

__excess_g_res(model,p,T,z,b,c) = excess_g_res(model,p,T,z)
function __excess_g_res(model::WilsonModel,p,T,z,b,c)
    V = diagvalues(b) .- c
    return excess_g_res_wilson(model,p,T,z,V)
end

function mixing_rule(model::PRModel,V,T,z,mixing_model::gErRuleModel,α,a,b,c)
    n = sum(z)
    #x = z./n
    invn = (one(n)/n)
    invn2 = invn^2
    gᴱᵣ = __excess_g_res(mixing_model.activity,1e5,T,z,b,c)
    b̄ = zero(gᴱᵣ)
    res = zero(T + first(z))
    for i in 1:length(model)
        zi,bi = z[i],b[i,i]
        zi2 = zi^2
        b̄ += bi*zi2
        res += zi*a[i,i]*α[i]/bi
        for j in 1:(i-1)
            zij = zi*z[j]
            b̄ += 2*b[i,j]*zij
        end
    end
    b̄ = b̄*invn2
    Λ = infinite_pressure_gibbs_correction(model,z)
    ā = (res + gᴱᵣ/Λ)*b̄*invn
    c̄ = dot(z,c)*invn
    return ā,b̄,c̄
end