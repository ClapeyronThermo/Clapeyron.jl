abstract type TVTPRTranslationModel <: TranslationModel end

struct TVTPRTranslationParam <: EoSParam
    Vc::SingleParam{Float64}
end

@newmodelsimple TVTPRTranslation TVTPRTranslationModel TVTPRTranslationParam


"""
    TVTPRTranslation <: TVTPRTranslationModel

    TVTPRTranslation(components;
    userlocations = String[],
    verbose::Bool=false)

## Model Parameters

- `Vc`: Single Parameter (`Float64`) - Critical Volume `[m³·mol⁻¹]`

## Description

Temperature-dependent VTPR Translation model for cubics:
```
V = V₀ + mixing_rule(cᵢ)

cᵢ = βᵢ*(Zc(PR) - Zcᵢ)*RTcᵢ/Pcᵢ
Zcᵢ = Pcᵢ*Vcᵢ/(RTcᵢ)
Zc(PR) = 0.30740130869870386

βᵢ = 0.35/(0.35+(ηᵢ*abs(Tr-αᵢ(T)))^γᵢ)
η = -74.458Zcᵢ + 26.966
γᵢ = 246.78Zcᵢ^2 - 107.21Zcᵢ + 12.67

```

## Model Construction Examples
```
# Using the default database
translation = TVTPRTranslation("water") #single input
translation = TVTPRTranslation(["water","ethanol"]) #multiple components

# Using user-provided parameters

# Passing files or folders
translation = TVTPRTranslation(["neon","hydrogen"]; userlocations = ["path/to/my/db","critical/Vc.csv"])

# Passing parameters directly
translation = TVTPRTranslation(["neon","hydrogen"];userlocations = (;Vc = [4.25e-5, 6.43e-5]))
```

## References

1. Ahlers, J., & Gmehling, J. (2001). Development of an universal group contribution equation of state. Fluid Phase Equilibria, 191(1–2), 177–188. [doi:10.1016/s0378-3812(01)00626-4](https://doi.org/10.1016/s0378-3812(01)00626-4) 
"""
TVTPRTranslation

export TVTPRTranslation

default_locations(::Type{TVTPRTranslation}) = critical_data()
default_references(::Type{TVTPRTranslation}) = ["10.1016/S0378-3812(01)00626-4"]

function translation(model::CubicModel,V,T,z,translation_model::TVTPRTranslation)
    res = zeros(eltype(V+T+first(z)),length(z))
    Tc = model.params.Tc.values
    Pc = model.params.Pc.values
    Vc = translation_model.params.Vc.values
    α = α_function(model,V,T,z,model.alpha)
    c = similar(z,Base.promote_eltype(model,T))
    #Zc_PR = Clapeyron.pure_cubic_zc(PR("water"))
    Zc_PR = 0.30740130869870386
    for i ∈ @comps
        Tci = Tc[i]
        Pci = Pc[i]
        RT = Tci*Rgas(model)
        Zc = Pci*Vc[i]/RT
        Tr = T/Tci
        η = -74.458*Zc+26.966
        γ = 246.78*Zc^2-107.21*Zc+12.67
        β = 0.35/(0.35+(η*abs(Tr-α[i]))^γ)
        cc = (Zc_PR - Zc)*RT/Pci
        c[i] = cc*β
    end
    return c
end

"""
    TVTPR(components;
    idealmodel = BasicIdeal,
    userlocations = String[],
    group_userlocations = String[]
    ideal_userlocations = String[],
    alpha_userlocations = String[],
    mixing_userlocations = String[],
    activity_userlocations = String[],
    translation_userlocations = String[],
    reference_state = nothing,
    verbose = false)

Volume-translated Peng Robinson equation of state with temperature-dependent translation. It uses the following models:
- Translation Model: [`TVTPRTranslation`](@ref)
- Alpha Model: [`TwuAlpha`](@ref)
- Mixing Rule Model: [`VTPRRule`](@ref) with [`VTPRUNIFAC`](@ref) activity

## References
1. Ahlers, J., & Gmehling, J. (2001). Development of an universal group contribution equation of state. Fluid Phase Equilibria, 191(1–2), 177–188. [doi:10.1016/s0378-3812(01)00626-4](https://doi.org/10.1016/s0378-3812(01)00626-4)
"""
function TVTPR(components;
    idealmodel = BasicIdeal,
    alpha = TwuAlpha, #here just for compatibility with the notebooks.
    translation = TVTPRTranslation,
    userlocations = String[], 
    group_userlocations = String[],
    ideal_userlocations = String[],
    alpha_userlocations = String[],
    mixing_userlocations = String[],
    activity_userlocations = String[],
    translation_userlocations = String[],
    reference_state = nothing,
    verbose = false)

    activity = VTPRUNIFAC(components,
            userlocations = activity_userlocations,
            group_userlocations = group_userlocations,
            verbose = verbose)

    _components = activity.groups.components #extract pure component list

    mixing = VTPRRule

    return PR(_components;
    idealmodel = idealmodel,
    alpha = alpha,
    mixing = mixing,
    activity = activity,
    translation = translation,
    userlocations = userlocations,
    ideal_userlocations = ideal_userlocations,
    alpha_userlocations = alpha_userlocations,
    mixing_userlocations = mixing_userlocations,
    translation_userlocations = translation_userlocations,
    reference_state = reference_state,
    verbose = verbose)
end

export TVTPR