"""
    EOS_LNG::MultiFluid
    EOS_LNG(components::Vector{String};Rgas = 8.314472,verbose = false)

## input Parameters

None

## Description

EOS-LNG: A Fundamental Equation of State for the Calculation of Thermodynamic Properties of Liquefied Natural Gases. valid for 21 compounds (`Clapeyron.GERG2008_names`). the EoS has new binary-specific parameters for methane + n-butane, methane + isobutane, methane + n-pentane, and methane + isopentane.

It uses the same functional form as [`GERG2008`](@ref).

## References

1. Thol, M., Richter, M., May, E. F., Lemmon, E. W., & Span, R. (2019). EOS-LNG: A fundamental equation of state for the calculation of thermodynamic properties of liquefied natural gases. Journal of Physical and Chemical Reference Data, 48(3), 033102. [doi:10.1063/1.5093800](https://doi.org/10.1063/1.5093800)
2. Kunz, O., & Wagner, W. (2012). The GERG-2008 wide-range equation of state for natural gases and other mixtures: An expansion of GERG-2004. Journal of Chemical and Engineering Data, 57(11), 3032–3091. [doi:10.1021/je300655b](https://doi.org/10.1021/je300655b)
"""
function EOS_LNG(components;verbose = false,Rgas = 8.314472)
    return MultiFluid(components;
    mixing = AsymmetricMixing,
    departure = EmpiricDeparture,
    pure_userlocations = String["@REMOVEDEFAULTS","@DB/Empiric/GERG2008/pures"],
    mixing_userlocations  = String["@REMOVEDEFAULTS","@DB/Empiric/EOS_LNG/mixing/EOS_LNG_mixing_unlike.csv"],
    departure_userlocations = String["@REMOVEDEFAULTS","@DB/Empiric/EOS_LNG/departure/EOS_LNG_departure_unlike.csv"],
    coolprop_userlocations = false,
    Rgas = Rgas,
    verbose = verbose)
end
export EOS_LNG
