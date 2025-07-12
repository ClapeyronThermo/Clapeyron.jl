abstract type PatelTejaModel <: ABCCubicModel end

struct PatelTejaParam <: EoSParam
    a::PairParam{Float64}
    b::PairParam{Float64}
    c::SingleParam{Float64}
    Tc::SingleParam{Float64}
    Pc::SingleParam{Float64}
    Vc::SingleParam{Float64}
    Mw::SingleParam{Float64}
end

function transform_params(::Type{PatelTejaParam},params,components)
    n = length(components)
    transform_params(ABCubicParam,params,components)
    Tc = params["Tc"]
    Pc = params["Pc"]
    c = get!(params,"c") do
        SingleParam("c",components,zeros(Base.promote_eltype(Pc,Tc),n),fill(true,n))
    end
    Vc = params["Vc"]


    Vc = get!(params,"Vc") do
        SingleParam("Vc",components,zeros(Base.promote_eltype(Pc,Tc),n),fill(true,n))
    end
    ω = get(params,"acentricfactor",nothing)
    if any(Vc.ismissingvalues)
        isnothing(ω) && throw(MissingException("PatelTeja: cannot estimate Vc: missing acentricfactor parameter."))
        
        for i in 1:n
            if Vc.ismissingvalues[i]
                ζc = evalpoly(ω[i],(0.329032,-0.076799,0.0211947))
                Vc[i] = ζc * R̄ * Tc[i] / Pc[i]
            end
        end
    end

    return params
end

struct PatelTeja{T <: IdealModel,α,c,γ} <:PatelTejaModel
    components::Array{String,1}
    alpha::α
    mixing::γ
    translation::c
    params::PatelTejaParam
    idealmodel::T
    references::Array{String,1}
end

"""
    PatelTeja(components;
    idealmodel = BasicIdeal,
    alpha = NoAlpha,
    mixing = vdW1fRule,
    activity = nothing,
    translation = PatelTejaTranslation,
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

Patel-Teja Equation of state.
```
P = RT/(v-b) + a•α(T)/((v - Δ₁b)*(v - Δ₂b))
aᵢᵢ = Ωaᵢ(R²Tcᵢ²/Pcᵢ)
bᵢᵢ = Ωbᵢ(R²Tcᵢ/Pcᵢ)
cᵢ = Ωcᵢ(R²Tcᵢ/Pcᵢ)
Zcᵢ = Pcᵢ*Vcᵢ/(R*Tcᵢ)
Ωaᵢ = 3Zcᵢ² + 3(1 - 2Zcᵢ)Ωbᵢ + Ωbᵢ² + 1 - 3Zcᵢ
0 = -Zcᵢ³ + (3Zcᵢ²)*Ωbᵢ + (2 - 3Zcᵢ)*Ωbᵢ² + Ωbᵢ³
Ωcᵢ = 1 - 3Zcᵢ

γ = ∑cᵢxᵢ/∑bᵢxᵢ
δ = 1 + 6γ + γ²
ϵ = 1 + γ

Δ₁ = -(ϵ + √δ)/2
Δ₂ = -(ϵ - √δ)/2
```
if `Vc` is not known, `Zc` can be estimated, using the acentric factor:

```
Zc = 0.329032 - 0.076799ω + 0.0211947ω²
```

## Model Construction Examples
```julia
# Using the default database
model = PatelTeja("water") #single input
model = PatelTeja(["water","ethanol"]) #multiple components
model = PatelTeja(["water","ethanol"], idealmodel = ReidIdeal) #modifying ideal model
model = PatelTeja(["water","ethanol"],alpha = TwuAlpha) #modifying alpha function
model = PatelTeja(["water","ethanol"],translation = RackettTranslation) #modifying translation
model = PatelTeja(["water","ethanol"],mixing = KayRule) #using another mixing rule
model = PatelTeja(["water","ethanol"],mixing = WSRule, activity = NRTL) #using advanced EoS+gᴱ mixing rule

# Passing a prebuilt model

my_alpha = PR78Alpha(["ethane","butane"],userlocations = Dict(:acentricfactor => [0.1,0.2]))
model = PatelTeja(["ethane","butane"],alpha = my_alpha)

# User-provided parameters, passing files or folders
model = PatelTeja(["neon","hydrogen"]; userlocations = ["path/to/my/db","cubic/my_k_values.csv"])

# User-provided parameters, passing parameters directly

model = PatelTeja(["neon","hydrogen"];
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

1. Patel, N. C., & Teja, A. S. (1982). A new cubic equation of state for fluids and fluid mixtures. Chemical Engineering Science, 37(3), 463–473. [doi:10.1016/0009-2509(82)80099-7](https://doi.org/10.1016/0009-2509(82)80099-7)

"""
PatelTeja

export PatelTeja
function PatelTeja(components;
    idealmodel = BasicIdeal,
    alpha = PatelTejaAlpha,
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
    model = CubicModel(PatelTeja,params,formatted_components;
                        idealmodel,alpha,mixing,activity,translation,
                        userlocations,ideal_userlocations,alpha_userlocations,activity_userlocations,mixing_userlocations,translation_userlocations,
                        reference_state, verbose)

    k = get(params,"k",nothing)
    l = get(params,"l",nothing)
    recombine_cubic!(model,k,l)
    set_reference_state!(model,reference_state;verbose)
    return model
end

default_references(::Type{PatelTeja}) = ["10.1016/0009-2509(82)80099-7"]

function ab_premixing(model::PatelTejaModel,mixing::MixingRule,k,l)
    _Tc = model.params.Tc
    _pc = model.params.Pc
    _Vc = model.params.Vc
    a = model.params.a
    b = model.params.b
    
    for i in 1:length(model)
        pci,Tci,Vci = _pc[i],_Tc[i],_Vc[i]
        Zc = pci * Vci / (R̄ * Tci)        
        poly = (-Zc^3,3Zc^2,2-3*Zc,1.0)
        _,Ωb1,Ωb2,Ωb3 = Solvers.real_roots3(poly)
        Ωb = max(Ωb1,Ωb2,Ωb3)
        Ωa = 3*Zc^2 + 3*(1 - 2*Zc)*Ωb + Ωb^2 + 1 - 3*Zc
        a[i] = Ωa*R̄^2*Tci^2/pci
        b[i] = Ωb*R̄*Tci/pci
    end
    epsilon_LorentzBerthelot!(a,k)
    sigma_LorentzBerthelot!(b,l)
    return a,b
end

function c_premixing(model::PatelTejaModel)
    _Tc = model.params.Tc
    _pc = model.params.Pc
    _Vc = model.params.Vc
    c = model.params.c
    _Zc = _pc .* _Vc ./ (R̄ .* _Tc)
    Ωc = @. 1 - 3*_Zc
    c .= Ωc .* R̄ .*_Tc ./ _pc
    return c
end

function cubic_Δ(model::PatelTejaModel,z)
    b = diagvalues(model.params.b.values)
    c = diagvalues(model.params.c.values)
    b̄ = dot(b,z)
    c̄ = dot(c,z)
    γ = c̄/b̄
    δ = sqrt(evalpoly(γ,(1,6,1)))
    ϵ = 1 + γ
    return (-0.5*(ϵ + δ), -0.5*(ϵ - δ))
end
