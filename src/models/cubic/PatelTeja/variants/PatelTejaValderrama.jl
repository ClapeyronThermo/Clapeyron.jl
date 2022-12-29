abstract type PTVModel <: PatelTejaModel end

const PTVParam = ABCCubicParam

struct PTV{T <: IdealModel,α,c,γ} <:PTVModel
    components::Array{String,1}
    alpha::α
    mixing::γ
    translation::c
    params::PTVParam
    idealmodel::T
    references::Array{String,1}
end

@registermodel PTV
"""
    PTV(components::Vector{String};
    idealmodel=BasicIdeal,
    alpha = NoAlpha,
    mixing = vdW1fRule,
    activity=nothing,
    translation=PTVTranslation,
    userlocations=String[],
    ideal_userlocations=String[],
    alpha_userlocations = String[],
    mixing_userlocations = String[],
    activity_userlocations = String[],
    translation_userlocations = String[],
    verbose=false)

## Input parameters
- `Tc`: Single Parameter (`Float64`) - Critical Temperature `[K]`
- `Pc`: Single Parameter (`Float64`) - Critical Pressure `[Pa]`
- `Vc`: Single Parameter (`Float64`) - Critical Volume `[m3/mol]`
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g/mol]`
- `k`: Pair Parameter (`Float64`) (optional)
- `l`: Pair Parameter (`Float64`) (optional)

## Model Parameters
- `Tc`: Single Parameter (`Float64`) - Critical Temperature `[K]`
- `Pc`: Single Parameter (`Float64`) - Critical Pressure `[Pa]`
- `Vc`: Single Parameter (`Float64`) - Critical Volume `[m3/mol]`
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g/mol]`
- `a`: Pair Parameter (`Float64`)
- `b`: Pair Parameter (`Float64`)
- `c`: Pair Parameter (`Float64`)

## Description

Patel-Teja-Valderrama Equation of state.

```
P = RT/(v-b) + a•α(T)/((v - Δ₁b)*(v - Δ₂b))
aᵢᵢ = Ωaᵢ(R²Tcᵢ²/Pcᵢ)
bᵢᵢ = Ωbᵢ(R²Tcᵢ/Pcᵢ)
cᵢ = Ωcᵢ(R²Tcᵢ/Pcᵢ)
Zcᵢ =  Pcᵢ*Vcᵢ/(R*Tcᵢ)
Ωaᵢ = 0.66121 - 0.76105Zcᵢ
Ωbᵢ = 0.02207 + 0.20868Zcᵢ
Ωcᵢ = 0.57765 - 1.87080Zcᵢ

γ = ∑cᵢxᵢ/∑bᵢxᵢ
δ = 1 + 6γ + γ²
ϵ = 1 + γ

Δ₁ =  -(ϵ + √δ)/2
Δ₂ =  -(ϵ - √δ)/2
```

## References

1. Valderrama, J. O. (1990). A generalized Patel-Teja equation of state for polar and nonpolar fluids and their mixtures. Journal of Chemical Engineering of Japan, 23(1), 87–91. [doi:10.1252/jcej.23.87](https://doi.org/10.1252/jcej.23.87)

"""
PTV

export PTV
function PTV(components::Vector{String}; idealmodel=BasicIdeal,
    alpha = PTVAlpha,
    mixing = vdW1fRule,
    activity=nothing,
    translation=NoTranslation,
    userlocations=String[],
    ideal_userlocations=String[],
    alpha_userlocations = String[],
    mixing_userlocations = String[],
    activity_userlocations = String[],
    translation_userlocations = String[],
     verbose=false)
    params = getparams(components, ["properties/critical.csv", "properties/molarmass.csv","SAFT/PCSAFT/PCSAFT_unlike.csv"]; userlocations=userlocations, verbose=verbose)
    k  = get(params,"k",nothing)
    l = get(params,"l",nothing)
    pc = params["Pc"]
    Vc = params["Vc"]
    Mw = params["Mw"]
    Tc = params["Tc"]
    init_mixing = init_model(mixing,components,activity,mixing_userlocations,activity_userlocations,verbose)
    a = PairParam("a",components,zeros(length(components)))
    b = PairParam("b",components,zeros(length(components)))
    c = PairParam("c",components,zeros(length(components)))
    init_idealmodel = init_model(idealmodel,components,ideal_userlocations,verbose)
    init_alpha = init_model(alpha,components,alpha_userlocations,verbose)
    init_translation = init_model(translation,components,translation_userlocations,verbose)
    packagedparams = PTVParam(a,b,c,Tc,pc,Vc,Mw)
    references = String["10.1252/jcej.23.87"]
    model = PTV(components,init_alpha,init_mixing,init_translation,packagedparams,init_idealmodel,references)
    recombine_cubic!(model,k,l)
    return model
end

function ab_premixing(model::PTVModel,mixing::MixingRule,k,l)
    _Tc = model.params.Tc
    _pc = model.params.Pc
    _Vc = model.params.Vc
    a = model.params.a
    b = model.params.b
    _Zc = @. _pc.*_Vc./(R̄*_Tc)
    Ωa = @. 0.66121-0.76105*_Zc
    Ωb = @. 0.02207+0.20868*_Zc
    diagvalues(a) .= @. Ωa*R̄^2*_Tc^2/_pc
    diagvalues(b) .= @. Ωb*R̄*_Tc/_pc
    epsilon_LorentzBerthelot!(a,k)
    sigma_LorentzBerthelot!(b,l)
    return a,b
end

function c_premixing(model::PTVModel)
    _Tc = model.params.Tc
    _pc = model.params.Pc
    _Vc = model.params.Vc
    c = model.params.c
    _Zc = @. _pc.*_Vc./(R̄*_Tc)

    Ωc = @. 0.57765-1.87080*_Zc
    diagvalues(c) .= Ωc .* R̄ .*_Tc ./ _pc
    c = sigma_LorentzBerthelot!(c)
    return c
end
#=
 (-B2-2(B2+B)+A)
 (-B2-2B2-2B+A)
 (-3B2-2B+A)
=#