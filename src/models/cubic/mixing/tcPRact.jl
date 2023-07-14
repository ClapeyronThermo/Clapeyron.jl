abstract type tcPRactRuleModel <: ActivityMixingRule end

struct tcPRactRule{γ} <: tcPRactRuleModel
    components::Array{String,1}
    activity::γ
    references::Array{String,1}
end

"""
    tcPRactRule{γ} <: tcPRactRuleModel

    tcPRactRule(components::Vector{String};
    activity = Wilson,
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
bᵢⱼ = (bᵢ^(2/3) + bⱼ^(2/3))^(3/2)/2
b̄ = ∑bᵢⱼxᵢxⱼ
c̄ = ∑cᵢxᵢ
ā = b̄RT*(∑xᵢaᵢ/bᵢ + gᴱᵣ/Λ
Λ = 1/(r₂ - r₁) * log((1 - r₂)/(1 - r₁))

if the model is Peng-Robinson:
    q = 0.53
if the model is Redlich-Kwong:
    q = 0.593
```

to use different values for `q`, overload `Clapeyron.tcPRactq(::CubicModel,::tcPRactModel) = q`


## References
1. Michelsen, M. L. (1990). A modified Huron-Vidal mixing rule for cubic equations of state. Fluid Phase Equilibria, 60(1–2), 213–219. [doi:10.1016/0378-3812(90)85053-d](https://doi.org/10.1016/0378-3812(90)85053-d)

"""
tcPRactRule

export tcPRactRule
function tcPRactRule(components::Vector{String}; activity = tcPRWilson, userlocations=String[],activity_userlocations=String[], verbose::Bool=false)
    _activity = init_model(activity,components,activity_userlocations,verbose)
    references = ["10.1016/0378-3812(90)85053-D"]
    model = tcPRactRule(components, _activity,references)
    return model
end

function ab_premixing(model::PRModel,mixing::tcPRactRuleModel,k = nothing, l = nothing)
    Ωa, Ωb = ab_consts(model)
    _Tc = model.params.Tc
    _pc = model.params.Pc
    a = model.params.a
    b = model.params.b
    diagvalues(a) .= @. Ωa*R̄^2*_Tc^2/_pc
    diagvalues(b) .= @. Ωb*R̄*_Tc/_pc
    epsilon_LorentzBerthelot!(a,k)
    tcPRact_mix(bi,bj,lij) = mix_powmean(bi,bj,lij,2/3)
    kij_mix!(tcPRact_mix,b,l)
    return a,b
end


__excess_g_res(model,p,T,z,b,c) = excess_g_res(model,p,T,z)
function __excess_g_res(model::WilsonModel,p,T,z,b,c)
    V = diagvalues(b) .- c
    return excess_g_res_wilson(model,p,T,z,V)
end

function mixing_rule(model::PRModel,V,T,z,mixing_model::tcPRactRuleModel,α,a,b,c)
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