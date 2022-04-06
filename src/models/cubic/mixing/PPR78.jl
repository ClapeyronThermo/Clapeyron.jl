abstract type PPR78RuleModel <: MixingRule end


struct PPR78Param <: EoSParam
    A::PairParam{Float64}
    B::PairParam{Float64}
end

struct PPR78Rule <: PPR78RuleModel
    groups::GroupParam
    components::Vector{String}
    params::PPR78Param
    references::Vector{String}
end
@registermodel HVRule

"""
    HVRule{γ} <: HVRuleModel
    
    HVRule(components;
    userlocations::Vector{String}=String[],
    verbose::Bool=false)

## Input Parameters

- `A`: Pair Parameter (`Float64`) - Fitted Parameter `[K]`
- `B`: Pair Parameter (`Float64`) - Fitted Parameter `[K]`

## Description

PPR78 Mixing Rule 
```
aᵢⱼ = √(aᵢaⱼ)
bᵢⱼ = (bᵢ +bⱼ)/2
b̄ = ∑bᵢⱼxᵢxⱼ
c̄ = ∑cᵢxᵢ
ā = b̄(∑[xᵢaᵢᵢαᵢ/(bᵢᵢ)] - ∑xᵢxⱼbᵢbⱼEᵢⱼ/2b̄)
Eᵢⱼ = ∑(z̄ᵢₖ - z̄ⱼₖ)(z̄ᵢₗ - z̄ⱼₗ) × Aₖₗ × (298.15/T)^(Aₖₗ/Bₖₗ - 1)
```

## References
1. Jaubert, J.-N., Privat, R., & Mutelet, F. (2010). Predicting the phase equilibria of synthetic petroleum fluids with the PPR78 approach. AIChE Journal. American Institute of Chemical Engineers, 56(12), 3225–3235. doi:10.1002/aic.12232
2. Jaubert, J.-N., Qian, J.-W., Lasala, S., & Privat, R. (2022). The impressive impact of including enthalpy and heat capacity of mixing data when parameterising equations of state. Application to the development of the E-PPR78 (Enhanced-Predictive-Peng-Robinson-78) model. Fluid Phase Equilibria, (113456), 113456. doi:10.1016/j.fluid.2022.113456

"""
HVRule

export HVRule
function HVRule(components::Vector{String}; activity = Wilson, userlocations::Vector{String}=String[],activity_userlocations::Vector{String}=String[], verbose::Bool=false)
    init_activity = activity(components;userlocations = activity_userlocations,verbose)
    
    references = ["10.1002/aic.12232","10.1016/j.fluid.2022.113456"]
    model = HVRule(components, init_activity,references)
    return model
end


function mixing_rule(model::CubicModel,V,T,z,mixing_model::HVRuleModel,α,a,b,c)
    n = sum(z)
    invn = 1/n
    invn2 = invn*invn
    b̄ = dot(z,Symmetric(b),z) * invn2
    c̄ = dot(z,c)/n
    gᴱ = zero(T+first(z))
    A = mixing_model.params.A.values
    B = mixing_model.params.B.values
    groups = mixing_model.groups
    gc = groups
    for i ∈ @comps
        for j in 1:i-1
            for k in gc
            end
        end
    end
    ∑ab = sum(z[i]*a[i,i]*α[i]/b[i,i] for i ∈ @comps)*invn
    ā = b̄*(∑ab-gᴱ/_λ)
    return ā,b̄,c̄
end

