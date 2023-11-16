abstract type WSRuleModel <: ActivityMixingRule end

struct WSRule{γ} <: WSRuleModel
    components::Array{String,1}
    activity::γ
    references::Array{String,1}
end

"""
    WSRule{γ} <: WSRuleModel

    WSRule(components;
    activity = Wilson,
    userlocations=String[],
    activity_userlocations=String[],
    verbose::Bool=false)

## Input Parameters

None

## Input models

- `activity`: Activity Model

## Description

Wong-Sandler Mixing Rule.

```
aᵢⱼ = √(aᵢaⱼ)(1 - kᵢⱼ)
bᵢⱼ = (bᵢ + bⱼ)/2
c̄ = ∑cᵢxᵢ
B̄ = ΣxᵢxⱼB̄ᵢⱼ
B̄ᵢⱼ = (1 - kᵢⱼ)((bᵢ - aᵢ/RT) + (bⱼ - aⱼ/RT))/2
b̄  = B̄/(1 - gᴱ/λRT - Σxᵢaᵢαᵢ/bᵢRT)
ā = RT(b̄ - B̄)
for Redlich-Kwong:
    λ = log(2) (0.6931471805599453)
for Peng-Robinson:
    λ = 1/(2√(2))log((2+√(2))/(2-√(2))) (0.6232252401402305)
```

`λ` is a coefficient indicating the relation between `gᴱ` and `gᴱ(cubic)` at infinite pressure. see [1] for more information. it can be customized by defining `WS_λ(::WSRuleModel,::CubicModel)`

## Model Construction Examples
```
# Using the default database
mixing = WSRule(["water","carbon dioxide"]) #default: Wilson Activity Coefficient.
mixing = WSRule(["water","carbon dioxide"],activity = NRTL) #passing another Activity Coefficient Model.
mixing = WSRule([("ethane",["CH3" => 2]),("butane",["CH2" => 2,"CH3" => 2])],activity = UNIFAC) #passing a GC Activity Coefficient Model.

# Passing a prebuilt model

act_model = NRTL(["water","ethanol"],userlocations = (a = [0.0 3.458; -0.801 0.0],b = [0.0 -586.1; 246.2 0.0], c = [0.0 0.3; 0.3 0.0]))
mixing = WSRule(["water","ethanol"],activity = act_model)

# Using user-provided parameters

# Passing files or folders
mixing = WSRule(["water","ethanol"]; activity = NRTL, activity_userlocations = ["path/to/my/db","nrtl_ge.csv"])

# Passing parameters directly
mixing = WSRule(["water","ethanol"];
                activity = NRTL,
                userlocations = (a = [0.0 3.458; -0.801 0.0],
                    b = [0.0 -586.1; 246.2 0.0],
                    c = [0.0 0.3; 0.3 0.0])
                )
```


## References
1. Wong, D. S. H., & Sandler, S. I. (1992). A theoretically correct mixing rule for cubic equations of state. AIChE journal. American Institute of Chemical Engineers, 38(5), 671–680. [doi:10.1002/aic.690380505](https://doi.org/10.1002/aic.690380505)
"""
WSRule

export WSRule
function WSRule(components; activity = Wilson, userlocations=String[],activity_userlocations=String[], verbose::Bool=false)
    _activity = init_mixing_act(activity,components,activity_userlocations,verbose)
    references = ["10.1002/aic.690380505"]
    model = WSRule(format_components(components), _activity,references)
    return model
end

WS_λ(::WSRuleModel,model::ABCubicModel,z) = infinite_pressure_gibbs_correction(model,z)

function mixing_rule(model::ABCubicModel,V,T,z,mixing_model::WSRuleModel,α,a,b,c)
    λ = WS_λ(mixing_model,model,z)
    n = sum(z)
    invn = (one(n)/n)
    RT⁻¹ = 1/(R̄*T)
    B̄ = zero(T+V+first(z))
    Σab = B̄
    for i in @comps
        zi = z[i]
        αi = α[i]
        _ai = a[i,i]
        ai = _ai*αi
        bi = b[i,i]
        B̄ += zi*zi*(bi-ai*RT⁻¹)
        Σab += zi*ai/bi
        for j in 1:(i-1)
            αj = α[j]
            bj= b[j,j]
            _aj = a[j,j]
            _1mkij = a[i,j]^2/(_ai*_aj) #1 - kij
            aj = a[j,j]*αj
            B̄ += _1mkij*zi*z[j]*((bj-aj*RT⁻¹)+(bi-ai*RT⁻¹))
        end
    end
    Σab = Σab*invn
    B̄ = B̄*invn*invn
    Aᴱ = excess_gibbs_free_energy(mixing_model.activity,1e5,T,z)*invn
    b̄  = B̄/(1 + (Aᴱ/λ - Σab)*RT⁻¹)
    ā = b̄*(Σab-Aᴱ/λ)
    c̄ = dot(z,c)*invn
    return ā,b̄,c̄
end

function cubic_get_k(model::CubicModel,mixing::WSRuleModel,params)
    return get_k_geomean(params.a.values)
end

function cubic_get_l(model::CubicModel,mixing::WSRuleModel,params)
    return get_k_mean(params.b.values)
end