const RKPRParam = ABCCubicParam
abstract type RKPRModel <: ABCCubicModel end

struct RKPR{T <: IdealModel,α,c,M} <: RKPRModel
    components::Array{String,1}
    alpha::α
    mixing::M
    translation::c
    params::RKPRParam
    idealmodel::T
    references::Array{String,1}
end

@registermodel RKPR
export RKPR

"""
    RKPR(components::Vector{String}; idealmodel=BasicIdeal,
    alpha = RKPRAlpha,
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

## Input models
- `idealmodel`: Ideal Model
- `alpha`: Alpha model
- `mixing`: Mixing model
- `activity`: Activity Model, used in the creation of the mixing model.
- `translation`: Translation Model

## Description
Redlich-Kwong-Peng-Robinson Equation of state.
```
P = RT/(v-b) + a•α(T)/((v + Δ₁b)*(v + Δ₂b))
Δ₁ = δ
Δ₂ = (1 - δ)/(1 + δ)
δ = ∑cᵢxᵢ

aᵢᵢ = Ωaᵢ(R²Tcᵢ²/Pcᵢ)
bᵢᵢ = Ωbᵢ(R²Tcᵢ/Pcᵢ)
Ωaᵢ = (3*yᵢ^2 + 3yᵢdᵢ + dᵢ^2 + dᵢ - 1)/(3yᵢ + dᵢ - 1)^2
Ωbᵢ = 1/(3yᵢ + dᵢ - 1)
dᵢ = (1 + cᵢ^2)/(1 + cᵢ)
yᵢ = 1 + (2(1 + cᵢ))^(1/3) + (4/(1 + cᵢ))^(1/3)
```
`cᵢ` is fitted to match:
```
if Zcᵢ[exp] > 0.29
    cᵢ =  √2 - 1
else
    Zcᵢ = 1.168Zcᵢ[exp]
    f(cᵢ) == 0
    f(cᵢ) = Zcᵢ - yᵢ/(3yᵢ + dᵢ - 1)
```
## References
1. Cismondi, M., & Mollerup, J. (2005). Development and application of a three-parameter RK–PR equation of state. Fluid Phase Equilibria, 232(1–2), 74–89. [doi:10.1016/j.fluid.2005.03.020](https://doi.org/10.1016/j.fluid.2005.03.020)
2. Tassin, N. G., Mascietti, V. A., & Cismondi, M. (2019). Phase behavior of multicomponent alkane mixtures and evaluation of predictive capacity for the PR and RKPR EoS’s. Fluid Phase Equilibria, 480, 53–65. [doi:10.1016/j.fluid.2018.10.005](https://doi.org/10.1016/j.fluid.2018.10.005)
"""
RKPR

function RKPR(components::Vector{String}; idealmodel=BasicIdeal,
    alpha = RKPRAlpha,
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
    n = length(components)
    a = PairParam("a",components,zeros(n))
    b = PairParam("b",components,zeros(n))
    c = PairParam("c",components,zeros(n))
    init_idealmodel = init_model(idealmodel,components,ideal_userlocations,verbose)
    init_alpha = init_model(alpha,components,alpha_userlocations,verbose)
    init_translation = init_model(translation,components,translation_userlocations,verbose)
    packagedparams = ABCCubicParam(a,b,c,Tc,pc,Vc,Mw)
    references = String["10.1016/j.fluid.2005.03.020","10.1016/j.fluid.2018.10.005"]
    model = RKPR(components,init_alpha,init_mixing,init_translation,packagedparams,init_idealmodel,references)
    recombine_cubic!(model,k,l)
    return model
end

function __rkpr_f0_δ(δ,Zc)
    δ2 = δ*δ
    d1 = (1 + δ2)/(1 + δ)
    y = 1 + cbrt(2*(1+δ)) + cbrt(4/(1+δ))
    return Zc - y/(3*y + d1 - 1)
end

function ab_premixing(model::RKPRModel,mixing::MixingRule,k  = nothing, l = nothing)
    _Tc = model.params.Tc
    _pc = model.params.Pc
    _Vc = model.params.Vc
    a = model.params.a
    b = model.params.b
    c = model.params.c
    prob = Roots.ZeroProblem(__rkpr_f0_δ,0.0)
    for i in @comps
        pci,Tci,Vci = _pc[i],_Tc[i],_Vc[i]
        Zci = pci * Vci / (R̄ * Tci)
        #Roots.find_zero(x -> Clapeyron.__rkpr_f0_δ(sqrt(2) - 1,1.168*x),0.29) 
        #0.2897160510687658
        if Zci >  0.2897160510687658
            δ = sqrt(2) - 1
        else
            Zci_eos = 1.168*Zci
            δ = Roots.solve(prob,Roots.Order0(),Zci_eos)
        end
        c[i] = δ
        d = (1 + δ*δ)/(1+δ)
        y = 1 + cbrt(2*(1+δ)) + cbrt(4/(1+δ))
        Ωa = (3*y*(y + d) + d*d + d - 1)/abs2(3*y + d - 1)
        Ωb = 1/(3*y + d - 1)
        a[i] = Ωa*R̄^2*Tci^2/pci
        b[i] = Ωb*R̄*Tci/pci
    end
    epsilon_LorentzBerthelot!(a,k)
    sigma_LorentzBerthelot!(b,l)
    return a,b
end

#premixing of c is done at ab_premixing level
c_premixing(model::RKPRModel) = model.params.c

function cubic_Δ(model::RKPRModel,z)
    c = diagvalues(model.params.c)
    z⁻¹ = sum(z)^-1
    Δ1 = zero(eltype(z))
    Δ2 = zero(Δ1)
    for i in @comps
        δi =  c[i]
        zi = z[i]
        Δ2 += zi*δi
        Δ1 += z[i]*((1 - δi)/(1 + δi))
    end
    return  -Δ2*z⁻¹, -Δ1*z⁻¹
end

function pure_cubic_zc(model::RKPRModel)
    δ = model.params.c.values[1]
    d = (1 + δ*δ)/(1+δ)
    y = 1 + cbrt(2*(1+δ)) + cbrt(4/(1+δ))
    Zc = y/(3y + d - 1)
    return Zc
end

crit_pure(model::RKPRModel) = crit_pure_tp(model)
function crit_pure_tp(model::RKPRModel)
    Tc = model.params.Tc.values[1]
    Pc = model.params.Pc.values[1]
    Zc = pure_cubic_zc(model) #PV = ZRT
    return (Tc,Pc,Zc*R̄*Tc/Pc)
end