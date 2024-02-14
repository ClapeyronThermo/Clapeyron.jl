"""
    GERG2008::MultiFluid

    GERG2008(components;
            Rgas = 8.314472,
            reference_state = nothing,
            verbose = false)

## input Parameters

None

## Description

The GERG-2008 Wide-Range Equation of State for Natural Gases and Other Mixtures. valid for 21 compounds (`Clapeyron.GERG2008_names`).
```

a = a⁰ + aʳ

a⁰ = ∑xᵢ(a⁰ᵢ(τᵢ,δᵢ) + ln(xᵢ))
δᵢ = ρ/ρcᵢ
τᵢ = Tcᵢ/T
a⁰ᵢ = ln(δᵢ) + R∗/R[n⁰ᵢ₋₁ + n⁰ᵢ₋₂τᵢ + n⁰ᵢ₋₃ln(τᵢ) + ∑n⁰ᵢ₋ₖln(abs(sinh(ϑ₀ᵢ₋ₖτᵢ))) + ∑n⁰ᵢ₋ₖln(cosh(ϑ₀ᵢ₋ₖτᵢ))]
R∗ = 8.314510
R = 8.314472

τ = Tᵣ/T
δ = ρ/ρᵣ
(1/ρᵣ) = ∑∑xᵢxⱼβᵥ₋ᵢⱼγᵥ₋ᵢⱼ[(xᵢ+xⱼ)/(xᵢβᵥ₋ᵢⱼ^2 + xⱼ)]•1/8(1/∛ρcᵢ + 1/∛ρcⱼ)^2
Tᵣ = ∑∑xᵢxⱼβₜ₋ᵢⱼγₜ₋ᵢⱼ[(xᵢ+xⱼ)/(xᵢβₜ₋ᵢⱼ^2 + xⱼ)]•√(TcᵢTcⱼ)
aʳ = ∑xᵢaᵣᵢ(τ,δ) + ∑∑xᵢxⱼFᵢⱼaʳᵢⱼ(τ,δ)
aʳᵢ = ∑nᵢ₋ₖδ^(dᵢ₋ₖ)τ^(tᵢ₋ₖ)  + ∑nᵢ₋ₖδ^(dᵢ₋ₖ)τ^(tᵢ₋ₖ)exp(-δ^cᵢ₋ₖ)
aʳᵢⱼ = ∑nᵢⱼ₋ₖδ^(dᵢⱼ₋ₖ)τ^(tᵢⱼ₋ₖ)  + ∑nᵢⱼ₋ₖδ^(dᵢⱼ₋ₖ)τ^(tᵢⱼ₋ₖ)exp(ηᵢⱼ₋ₖ(δ-εᵢⱼ₋ₖ)^2 + βᵢⱼ₋ₖ(δ-γᵢⱼ₋ₖ))
```

## References

1. Kunz, O., & Wagner, W. (2012). The GERG-2008 wide-range equation of state for natural gases and other mixtures: An expansion of GERG-2004. Journal of Chemical and Engineering Data, 57(11), 3032–3091. [doi:10.1021/je300655b](https://doi.org/10.1021/je300655b)
"""
function GERG2008(components;            
    Rgas = 8.314472,
    reference_state = nothing,
    verbose = false)

    return MultiFluid(components;
    mixing = AsymmetricMixing,
    departure = EmpiricDeparture,
    pure_userlocations = String["@REMOVEDEFAULTS","@DB/Empiric/GERG2008/pures"],
    mixing_userlocations  = String["@REMOVEDEFAULTS","@DB/Empiric/GERG2008/mixing/GERG2008_mixing_unlike.csv"],
    departure_userlocations = String["@REMOVEDEFAULTS","@DB/Empiric/GERG2008/departure/GERG2008_departure_unlike.csv"],
    coolprop_userlocations = false,
    Rgas = Rgas,
    reference_state = reference_state,
    verbose = verbose)
end
export GERG2008


function ref_f()
    model = GERG2008("water",reference_state = ReferenceState(:nbp))
    T,v,_ = saturation_temperature(model,101325.0)
    h = Clapeyron.VT_enthalpy(model,v,T)
    s = Clapeyron.VT_entropy(model,v,T)
    return (;h,s)
end