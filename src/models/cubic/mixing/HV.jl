abstract type HVRuleModel <: ActivityMixingRule end

struct HVRule{γ} <: HVRuleModel
    components::Array{String,1}
    activity::γ
    references::Array{String,1}
end
@registermodel HVRule

"""
    HVRule{γ} <: HVRuleModel

    HVRule(components::Vector{String};
    activity = Wilson,
    userlocations::Vector{String}=String[],
    activity_userlocations::Vector{String}=String[],
    verbose::Bool=false)

## Input Parameters

None

## Input models

- `activity`: Activity Model

## Description

Huron-Vidal Mixing Rule
```
aᵢⱼ = √(aᵢaⱼ)(1-kᵢⱼ)
bᵢⱼ = (bᵢ +bⱼ)/2
b̄ = ∑bᵢⱼxᵢxⱼ
c̄ = ∑cᵢxᵢ
ā = b̄(∑[xᵢaᵢᵢαᵢ/(bᵢᵢ)] - gᴱ/λ)
for Redlich-Kwong:
    λ = log(2) (0.6931471805599453)
for Peng-Robinson:
    λ = 1/(2√(2))log((2+√(2))/(2-√(2))) (0.6232252401402305)
```
`λ` is a coefficient indicating the relation between `gᴱ` and `gᴱ(cubic)` at infinite pressure. see [1] for more information. it can be customized by defining `HV_λ(::HVRuleModel,::CubicModel,z)`

## References
1. Huron, M.-J., & Vidal, J. (1979). New mixing rules in simple equations of state for representing vapour-liquid equilibria of strongly non-ideal mixtures. Fluid Phase Equilibria, 3(4), 255–271. [doi:10.1016/0378-3812(79)80001-1](https://doi.org/10.1016/0378-3812(79)80001-1)

"""
HVRule

export HVRule
function HVRule(components::Vector{String}; activity = Wilson, userlocations::Vector{String}=String[],activity_userlocations::Vector{String}=String[], verbose::Bool=false)
    _activity = init_model(activity,components,activity_userlocations,verbose)
    references = ["10.1016/0378-3812(79)80001-1"]
    model = HVRule(components, _activity,references)
    return model
end

HV_λ(::HVRuleModel,model::ABCubicModel,z) = infinite_pressure_gibbs_correction(model,z)

function mixing_rule(model::ABCubicModel,V,T,z,mixing_model::HVRuleModel,α,a,b,c)
    n = sum(z)
    invn = 1/n
    invn2 = invn*invn
    b̄ = dot(z,Symmetric(b),z) * invn2
    c̄ = dot(z,c)/n
    gᴱ = excess_gibbs_free_energy(mixing_model.activity,1e5,T,z)*invn
    ∑ab = sum(z[i]*a[i,i]*α[i]/b[i,i] for i ∈ @comps)*invn
    _λ = HV_λ(mixing_model,model,z)
    ā = b̄*(∑ab-gᴱ/_λ)
    return ā,b̄,c̄
end

