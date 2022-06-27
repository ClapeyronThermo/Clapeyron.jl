abstract type VTPRRuleModel <: MixingRule end

struct VTPRRule{γ} <: VTPRRuleModel
    components::Array{String,1}
    activity::γ
    references::Array{String,1}
end

@registermodel VTPRRule

"""
    VTPRRule{γ} <: VTPRRuleModel
    
    VTPRRule(components::Vector{String};
    activity = UNIFAC,
    userlocations::Vector{String}=String[],
    activity_userlocations::Vector{String}=String[],
    verbose::Bool=false)

## Input Parameters

None

## Input models 

- `activity`: Activity Model

## Description

Mixing Rule used by the Volume-translated Peng-Robinson [`VTPR`](@ref) equation of state.
only works with activity models that define an excess residual gibbs energy function `Clapeyron.excess_g_res(model,P,T,z)` function (like [`UNIQUAC`](@ref) and [`UNIFAC`](@ref) models)

```
aᵢⱼ = √(aᵢaⱼ)(1-kᵢⱼ)
bᵢⱼ = ((bᵢ^(3/4) + bⱼ^(3/4))/2)^(4/3)
log(γʳ)ᵢ = lnγ_res(model.activity,V,T,z) 
gᴱᵣₑₛ = ∑RTlog(γʳ)ᵢxᵢ
b̄ = ∑bᵢⱼxᵢxⱼ
c̄ = ∑cᵢxᵢ
ā = b̄RT(∑[xᵢaᵢᵢαᵢ/(RTbᵢᵢ)] - gᴱᵣₑₛ/(0.53087RT))
```

## References
1. Ahlers, J., & Gmehling, J. (2001). Development of an universal group contribution equation of state. Fluid Phase Equilibria, 191(1–2), 177–188. [doi:10.1016/s0378-3812(01)00626-4](https://doi.org/10.1016/s0378-3812(01)00626-4)
"""
VTPRRule

export VTPRRule
function VTPRRule(components::Vector{String}; activity = UNIFAC, userlocations::Vector{String}=String[],activity_userlocations::Vector{String}=String[], verbose::Bool=false)
    init_activity = activity(components;userlocations = activity_userlocations,verbose)
    
    references = ["10.1016/S0378-3812(01)00626-4"]
    model = VTPRRule(components, init_activity,references)
    return model
end

function ab_premixing(::Type{<:PRModel},mixing::VTPRRule,Tc,pc,kij)
    Ωa, Ωb = ab_consts(PR)
    _Tc = Tc.values
    _pc = pc.values
    components = pc.components
    a = epsilon_LorentzBerthelot(SingleParam(pc, @. Ωa*R̄^2*_Tc^2/_pc),kij)
    bi = @. Ωb*R̄*_Tc/_pc
    vtpr_mix(bi,bj,kij) = mix_powmean(bi,bj,0,3/4)
    b = kij_mix(vtpr_mix,SingleParam("b (covolume)",components,bi))
    return a,b
end

function mixing_rule(model::PRModel,V,T,z,mixing_model::VTPRRuleModel,α,a,b,c)
    n = sum(z)
    invn = (one(n)/n)
    invn2 = invn^2
    g_E_res = excess_gibbs_free_energy(mixing_model.activity,1e5,T,z) / n
    b̄ = dot(z,Symmetric(b),z) * invn2
    c̄ = dot(z,c)/n
    Σab = invn*sum(z[i]*a[i,i]*α[i]/b[i,i]/(R̄*T) for i ∈ @comps)
    ā = b̄*R̄*T*(Σab-1/0.53087*(g_E_res/(R̄*T)))
    return ā,b̄,c̄
end