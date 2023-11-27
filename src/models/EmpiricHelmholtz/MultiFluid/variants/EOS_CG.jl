"""
    EOS_CG::MultiFluid
    EOS_CG(components::Vector{String};Rgas = 8.314472,verbose = false)

## Input Parameters

None

## Description

EOS-CG: : A Mixture Model for the Calculation of Thermodynamic Properties of CCS Mixtures

It uses the same functional form as [`GERG2008`](@ref).

## References

1. Neumann, T., Herrig, S., Bell, I. H., Beckmüller, R., Lemmon, E. W., Thol, M., & Span, R. (2023). EOS-CG-2021: A mixture model for the calculation of thermodynamic properties of CCS mixtures. International Journal of Thermophysics, 44(12). [doi:10.1007/s10765-023-03263-6](https://doi.org/10.1007/s10765-023-03263-6)
"""
function EOS_CG(components;verbose = false,Rgas = Rgas())
    return MultiFluid(components;
    mixing = AsymmetricMixing,
    departure = EmpiricDeparture,
    pure_userlocations = String["@REMOVEDEFAULTS","@DB/Empiric/EOS_CG/pures"],
    mixing_userlocations  = String["@REMOVEDEFAULTS","@DB/Empiric/EOS_CG/mixing/EOS_CG_mixing_unlike.csv"],
    departure_userlocations = String["@REMOVEDEFAULTS","@DB/Empiric/EOS_CG/departure/EOS_CG_departure_unlike.csv"],
    coolprop_userlocations = false,
    Rgas = Rgas,
    verbose = verbose)
end
export EOS_LNG