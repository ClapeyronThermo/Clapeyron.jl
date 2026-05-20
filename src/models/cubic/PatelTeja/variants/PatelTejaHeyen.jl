@newmodelsimple PTHAlpha TwuAlphaModel TwuAlphaParam
default_locations(::Type{PTHAlpha}) = critical_data()
default_references(::Type{PTHAlpha}) = ["10.1016/0378-3812(80)80003-3"]
default_ignore_missing_singleparams(::Type{PTHAlpha}) = ["Vc","L","M","N","acentricfactor","Tc","Pc"]

function transform_params(::Type{PTHAlpha},params,components)
    nc = length(components)
    M = get!(params,"M") do
        SingleParam("M",components,zeros(nc),fill(true,nc))
    end
    w = get!(params,"acentricfactor") do
        SingleParam("acentricfactor",components,zeros(nc),fill(true,nc))
    end

    N = get!(params,"N") do
        SingleParam("N",components,zeros(nc),fill(true,nc))
    end

    L = get!(params,"M") do
        SingleParam("M",components,zeros(nc),fill(true,nc))
    end

    for i in 1:nc
        if M.ismissingvalues[i]
            M[i] = 1
        end
    end

    for i in 1:nc
        if N.ismissingvalues[i] && L.ismissingvalues[i] && w.ismissingvalues[i]
            throw(error("PTHAlpha: cannot estimate M and L: missing acentricfactor parameter"))
        end
        
        if N.ismissingvalues[i] && L.ismissingvalues[i] && !(w.ismissingvalues[i])
            #=
            Heyen: alpha = exp(H1*(1 - Tr^H2))
            Twu:   alpha = Tr^(N*(M-1))*exp(L*(1 - Tr^M*N))
            
            we can convert from Heyen to Twu:
            L = H1
            N = H2
            M = 1 
            =#
            ω = w[i]
            βc = evalpoly(ω,(1.4563,1.26,−0.3928))
            γc = evalpoly(ω,(−0.2981,−1.9574,0.1789))
            H2 = γc/(βc - 1) + βc
            H1 = (βc - 1)/H2
            L[i] = H1
            N[i] = H2
        end
    end

    return params
end

struct PatelTejaHeyen{T <: IdealModel,α,c,γ} <:PatelTejaModel
    components::Array{String,1}
    alpha::α
    mixing::γ
    translation::c
    params::PatelTejaParam
    idealmodel::T
    references::Array{String,1}
end

"""
    PatelTejaHeyen(components;
    idealmodel = BasicIdeal,
    alpha = PTHAlpha,
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
- `Vc_fit`: Single Parameter (`Float64`) - Fitted Critical Volume `[m³·mol⁻¹]`
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g·mol⁻¹]`
- `k`: Pair Parameter (`Float64`) (optional)
- `l`: Pair Parameter (`Float64`) (optional)
- `acentricfactor`: Single Parameter (`Float64`) (optional) - acentric factor, used to fit the critical volume if no value is provided.


## Model Parameters
- `Tc`: Single Parameter (`Float64`) - Critical Temperature `[K]`
- `Pc`: Single Parameter (`Float64`) - Critical Pressure `[Pa]`
- `Vc_fit`: Single Parameter (`Float64`) - Fitted Critical Volume `[m³·mol⁻¹]`
- `acentricfactor`: Single Parameter (`Float64`)
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g·mol⁻¹]`
- `a`: Pair Parameter (`Float64`)
- `b`: Pair Parameter (`Float64`)
- `c`: Pair Parameter (`Float64`)

## Description

Patel-Teja-Heyen Equation of state.
```
P = RT/(v-b) + a•α(T)/((v - Δ₁b)*(v - Δ₂b))
aᵢᵢ = Ωaᵢ(R²Tcᵢ²/Pcᵢ)
bᵢᵢ = Ωbᵢ(R²Tcᵢ/Pcᵢ)
cᵢ = Ωcᵢ(R²Tcᵢ/Pcᵢ)
Zcᵢ = Pcᵢ*Vcᵢ/(R*Tcᵢ)
Ωaᵢ = 3Zcᵢ² + 3(1 - 2Zcᵢ)Ωbᵢ + Ωbᵢ² + 1 - 3Zcᵢ
Ωbᵢ: maximum real solution of 0 = -Zcᵢ³ + (3Zcᵢ²)*Ωbᵢ + (2 - 3Zcᵢ)*Ωbᵢ² + Ωbᵢ³
Ωcᵢ = 1 - 3Zcᵢ

γ = ∑cᵢxᵢ/∑bᵢxᵢ
δ = 1 + 6γ + γ²
ϵ = 1 + γ

Δ₁ = -(ϵ + √δ)/2
Δ₂ = -(ϵ - √δ)/2
```
if `Vc_fit` is not known, `Zc` can be estimated, using the acentric factor:

```
Zc = 0.3277 - 0.0537*ω - 0.0147*ω²
```

## Model Construction Examples
```julia
# Using the default database
model = PatelTejaHeyen("water") #single input
model = PatelTejaHeyen(["water","ethanol"]) #multiple components
model = PatelTejaHeyen(["water","ethanol"], idealmodel = ReidIdeal) #modifying ideal model
model = PatelTejaHeyen(["water","ethanol"],alpha = TwuAlpha) #modifying alpha function
model = PatelTejaHeyen(["water","ethanol"],translation = RackettTranslation) #modifying translation
model = PatelTejaHeyen(["water","ethanol"],mixing = KayRule) #using another mixing rule
model = PatelTejaHeyen(["water","ethanol"],mixing = WSRule, activity = NRTL) #using advanced EoS+gᴱ mixing rule

# Passing a prebuilt model

my_alpha = PR78Alpha(["ethane","butane"],userlocations = Dict(:acentricfactor => [0.1,0.2]))
model = PatelTejaHeyen(["ethane","butane"],alpha = my_alpha)

# User-provided parameters, passing files or folders
model = PatelTejaHeyen(["neon","hydrogen"]; userlocations = ["path/to/my/db","cubic/my_k_values.csv"])

# User-provided parameters, passing parameters directly

model = PatelTejaHeyen(["neon","hydrogen"];
        userlocations = (;Tc = [44.492,33.19],
                        Pc = [2679000, 1296400],
                        Mw = [20.17, 2.],
                        acentricfactor = [-0.03,-0.21]
                        k = [0. 0.18; 0.18 0.], #k,l can be ommited in single-component models.
                        l = [0. 0.01; 0.01 0.])
                    )
```

## References

1. Patel, N. C., & Teja, A. S. (1982). A new cubic equation of state for fluids and fluid mixtures. Chemical Engineering Science, 37(3), 463–473. [doi:10.1016/0009-2509(82)80099-7](https://doi.org/10.1016/0009-2509(82)80099-7)

"""
PatelTejaHeyen

export PatelTejaHeyen

function PatelTejaHeyen(components;
    idealmodel = BasicIdeal,
    alpha = PTHAlpha,
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
    params = getparams(formatted_components, ["properties/critical.csv", "properties/molarmass.csv","cubic/PatelTejaHeyen/PatelTeja_Vc_fit.csv"]; userlocations = userlocations, verbose = verbose,ignore_missing_singleparams = ["Vc_fit","acentricfactor"])
    model = CubicModel(PatelTejaHeyen,params,formatted_components;
                        idealmodel,alpha,mixing,activity,translation,
                        userlocations,ideal_userlocations,alpha_userlocations,activity_userlocations,mixing_userlocations,translation_userlocations,
                        reference_state, verbose)

    k = get(params,"k",nothing)
    l = get(params,"l",nothing)
    recombine_cubic!(model,k,l)
    set_reference_state!(model,reference_state;verbose)
    return model
end

function c_premixing(model::PatelTejaHeyen)
    _Tc = model.params.Tc
    _pc = model.params.Pc
    _Vc = model.params.Vc_fit
    ω = model.params.acentricfactor
    c = model.params.c

    for i in 1:length(model)
        Tc,Pc = _Tc[i],_pc[i]
        if _Vc.ismissingvalues[i]
            ω.ismissingvalues[i] && throw(MissingException("PatelTeja: cannot estimate Vc: missing acentricfactor parameter."))
            ζc = evalpoly(ω[i],(0.3277,-0.0537,-0.0147))
            _Vc[i] = ζc * R̄ * Tc / Pc
        end
        Zc = Pc*_Vc[i]/(R̄*Tc)
        Ωc = 1 - 3*Zc
        c[i] = Ωc*R̄*Tc/Pc
    end
    return c
end
