abstract type RKPRModel <: ABCubicModel end

struct RKPRParam <: EoSParam
    a::PairParam{Float64}
    b::PairParam{Float64}
    delta::SingleParam{Float64}
    Tc::SingleParam{Float64}
    Pc::SingleParam{Float64}
    Vc::SingleParam{Float64}
    Mw::SingleParam{Float64}
end

function transform_params(::Type{RKPRParam},params,components)
    n = length(components)
    transform_params(ABCubicParam,params,components)
    Tc = params["Tc"]
    Pc = params["Pc"]
    Vc = get(params,"Vc",1.0)
    delta = get!(params,"delta") do
        SingleParam("delta",components,zeros(Base.promote_eltype(Pc,Tc,Vc),n),fill(true,n))
    end
    return params
end

struct RKPR{T <: IdealModel,α,c,M} <: RKPRModel
    components::Array{String,1}
    alpha::α
    mixing::M
    translation::c
    params::RKPRParam
    idealmodel::T
    references::Array{String,1}
end

export RKPR

"""
    RKPR(components;
    idealmodel = BasicIdeal,
    alpha = RKPRAlpha,
    mixing = vdW1fRule,
    activity = nothing,
    translation = NoTranslation,
    userlocations = String[],
    ideal_userlocations = String[],
    alpha_userlocations = String[],
    mixing_userlocations = String[],
    activity_userlocations = String[],
    translation_userlocations = String[],
    reference_state = nothing,
    verbose = false)

## Input parameters
- `Tc`: Single Parameter (`Float64`) - Critical Temperature `[K]`
- `Pc`: Single Parameter (`Float64`) - Critical Pressure `[Pa]`
- `Vc`: Single Parameter (`Float64`) (optional) - Critical Volume `[m3/mol]`
- `delta`: Single Parameter (`Float64`) (optional) - constant parameter (no units)
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
- `delta`: Single Parameter (`Float64`) - constant parameter (no units)

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

If `c` is not provided, `cᵢ` is fitted to match:

```
if Zcᵢ[exp] > 0.29
    cᵢ = √2 - 1
else
    Zcᵢ = 1.168Zcᵢ[exp]
    f(cᵢ) == 0
    f(cᵢ) = Zcᵢ - yᵢ/(3yᵢ + dᵢ - 1)
```

## Model Construction Examples
```julia
# Using the default database
model = RKPR("water") #single input
model = RKPR(["water","ethanol"]) #multiple components
model = RKPR(["water","ethanol"], idealmodel = ReidIdeal) #modifying ideal model
model = RKPR(["water","ethanol"],alpha = TwuAlpha) #modifying alpha function
model = RKPR(["water","ethanol"],translation = RackettTranslation) #modifying translation
model = RKPR(["water","ethanol"],mixing = KayRule) #using another mixing rule
model = RKPR(["water","ethanol"],mixing = WSRule, activity = NRTL) #using advanced EoS+gᴱ mixing rule

# Passing a prebuilt model

my_alpha = PR78Alpha(["ethane","butane"],userlocations = Dict(:acentricfactor => [0.1,0.2]))
model = RKPR(["ethane","butane"],alpha = my_alpha)

# User-provided parameters, passing files or folders
model = RKPR(["neon","hydrogen"]; userlocations = ["path/to/my/db","cubic/my_k_values.csv"])

# User-provided parameters, passing parameters directly

model = RKPR(["neon","hydrogen"];
        userlocations = (;Tc = [44.492,33.19],
                        Pc = [2679000, 1296400],
                        Vc = [4.25e-5, 6.43e-5],
                        Mw = [20.17, 2.],
                        acentricfactor = [-0.03,-0.21]
                        k = [0. 0.18; 0.18 0.], #k,l can be ommited in single-component models.
                        l = [0. 0.01; 0.01 0.])
                    )
```
## References
1. Cismondi, M., & Mollerup, J. (2005). Development and application of a three-parameter RK–PR equation of state. Fluid Phase Equilibria, 232(1–2), 74–89. [doi:10.1016/j.fluid.2005.03.020](https://doi.org/10.1016/j.fluid.2005.03.020)
2. Tassin, N. G., Mascietti, V. A., & Cismondi, M. (2019). Phase behavior of multicomponent alkane mixtures and evaluation of predictive capacity for the PR and RKPR EoS’s. Fluid Phase Equilibria, 480, 53–65. [doi:10.1016/j.fluid.2018.10.005](https://doi.org/10.1016/j.fluid.2018.10.005)
"""
RKPR

function RKPR(components;
    idealmodel = BasicIdeal,
    alpha = RKPRAlpha,
    mixing = vdW1fRule,
    activity = nothing,
    translation = NoTranslation,
    userlocations = String[],
    ideal_userlocations = String[],
    alpha_userlocations = String[],
    mixing_userlocations = String[],
    activity_userlocations = String[],
    translation_userlocations = String[],
    reference_state = nothing,
    verbose = false)
    
    formatted_components = format_components(components)
    params = getparams(formatted_components, ["properties/critical.csv", "properties/molarmass.csv","SAFT/PCSAFT/PCSAFT_unlike.csv"]; userlocations = userlocations, verbose = verbose,ignore_missing_singleparams = ["Vc"])
    
    model = CubicModel(RKPR,params,formatted_components;
                        idealmodel,alpha,mixing,activity,translation,
                        userlocations,ideal_userlocations,alpha_userlocations,activity_userlocations,mixing_userlocations,translation_userlocations,
                        reference_state, verbose)

    k = get(params,"k",nothing)
    l = get(params,"l",nothing)
    recombine_cubic!(model,k,l)
    set_reference_state!(model,reference_state;verbose)
    return model
end

default_references(::Type{RKPR}) = ["10.1016/j.fluid.2005.03.020","10.1016/j.fluid.2018.10.005"]

function __rkpr_f0_δ(δ,Zc)
    δ2 = δ*δ
    d1 = (1 + δ2)/(1 + δ)
    y = 1 + cbrt(2*(1+δ)) + cbrt(4/(1+δ))
    return Zc - y/(3*y + d1 - 1)
end

function rkpr(δ)
    Δ1 = -δ
    Δ2 = -(1 - δ)/(1 + δ)
    return cubic_pure_zc(Δ1,Δ2)
end

function ab_premixing(model::RKPRModel,mixing::MixingRule,k, l)
    _Tc = model.params.Tc
    _pc = model.params.Pc
    _Vc = model.params.Vc
    a = model.params.a
    b = model.params.b
    delta = model.params.delta
    for i in @comps
        pci,Tci,Vci = _pc[i],_Tc[i],_Vc[i]
        δ = delta[i]
        Δ1 = -δ
        Δ2 = -(1 - δ)/(1 + δ)
        Ωa,Ωb = ab_consts(Δ1,Δ2)
        a[i] = Ωa*R̄^2*Tci^2/pci
        b[i] = Ωb*R̄*Tci/pci
    end
    epsilon_LorentzBerthelot!(a,k)
    sigma_LorentzBerthelot!(b,l)
    return a,b
end

function c_premixing(model::RKPRModel)
    _Tc = model.params.Tc
    _pc = model.params.Pc
    _Vc = model.params.Vc
    delta = model.params.delta
    prob = Roots.ZeroProblem(__rkpr_f0_δ,0.41*oneunit(eltype(delta))) #TODO: find a more stable way to solve this
    for i in @comps
        pci,Tci,Vci = _pc[i],_Tc[i],_Vc[i]
        if !_Vc.ismissingvalues[i] && delta.ismissingvalues[i]
            Zci = pci * Vci / (R̄ * Tci)
            #Roots.find_zero(x -> Clapeyron.__rkpr_f0_δ(sqrt(2) - 1,1.168*x),0.29)
            #0.2897160510687658
            if Zci > 0.2897160510687658
                δ = sqrt(2) - 1
            else
                Zci_eos = 1.168*Zci
                δ = Roots.solve(prob,Roots.Order0(),Zci_eos)
            end
            delta[i] = δ
        elseif !delta.ismissingvalues[i]
            
        else
            throw(MissingException("RKPR: Vc or delta needs to be specified."))
        end
    end
    return delta
end

function cubic_Δ(model::RKPRModel,z)
    δ = diagvalues(model.params.delta)
    z⁻¹ = sum(z)^-1
    Δ1 = zero(eltype(z))
    Δ2 = zero(Δ1)
    for i in @comps
        δi = δ[i]
        zi = z[i]
        Δ2 += zi*δi
        Δ1 += z[i]*((1 - δi)/(1 + δi))
    end
    return  -Δ2*z⁻¹, -Δ1*z⁻¹
end
