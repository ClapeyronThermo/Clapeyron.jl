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
    userlocations = String[],
    activity_userlocations = String[],
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

`λ` is a coefficient indicating the relation between `gᴱ` and `gᴱ(cubic)` at infinite pressure. See [1] for more information. It can be customized by defining `WS_λ(::WSRuleModel,::CubicModel)`.

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
function WSRule(components; activity = Wilson, userlocations = String[],activity_userlocations = String[], verbose::Bool=false)
    _activity = init_mixing_act(activity,components,activity_userlocations,verbose)
    references = ["10.1002/aic.690380505"]
    model = WSRule(format_components(components), _activity,references)
    return model
end

WS_λ(::WSRuleModel,model::DeltaCubicModel,T,z) = infinite_pressure_gibbs_correction(model,T,z)

function mixing_rule(model::DeltaCubicModel,V,T,z,mixing_model::WSRuleModel,α,a,b)
    λ = WS_λ(mixing_model,model,T,z)
    n = sum(z)
    nc = length(model)
    invn = (one(n)/n)
    RT⁻¹ = 1/(R̄*T)
    B̄ = zero(T+V+first(z))
    Σλab = B̄
    for i in 1:nc
        zi,αi,a0i,bi = z[i],α[i],a[i,i],b[i,i]
        ai = a0i*αi
        λi = WS_λ(mixing_model,model,T,FillArrays.OneElement(i,nc))
        Bi = (bi - ai*RT⁻¹)
        B̄ += zi*zi*(bi-ai*RT⁻¹)
        Σλab += λi*zi*ai/bi
        for j in 1:(i-1)
            zj,αj,a0j,bj = z[j],α[j],a[j,j],b[j,j]
            aj,aij = a0j*αj,a[i,j]
            k̄ij = aij*aij/(a0i*a0j) #1 - kij
            Bj = (bj - aj*RT⁻¹)
            B̄ += k̄ij*zi*zj*(Bi + Bj)
        end
    end
    Σλab = Σλab*invn
    B̄ = B̄*invn*invn
    Aᴱ = excess_gibbs_free_energy(mixing_model.activity,1e5,T,z)*invn
    b̄  = B̄/(1 + (Aᴱ - Σλab)/λ * RT⁻¹)
    ā = b̄*(Σλab-Aᴱ)/λ
    c̄ = translation2(model,V,T,z,model.translation,a,b,α)*invn
    return ā,b̄,c̄
end

function cubic_get_k(model::CubicModel,mixing::WSRuleModel,params)
    return get_k_geomean(params.a.values)
end

function cubic_get_l(model::CubicModel,mixing::WSRuleModel,params)
    return get_k_mean(params.b.values)
end
