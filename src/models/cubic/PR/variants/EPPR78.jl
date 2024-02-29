#=
groups:
CH3
CH2
CH
C
CH4
C2H6
CH aro
C aro
C fused aromatic rings
CH2 cyclic
CH cyclic~|~C cyclic
CO2
N2
H2S
SH
H2O
C2H4
CH2 alkenic~|~CH alkenic
C alkenic
CH cycloalkenic~|~C cycloalkenic
H2
C2F6
CF3
CF2
CF2 double bond~|~CF double bond
C2H4F2
C2H2F4
CO
He
Ar
SO2
O2
NO
COS
NH3
NO2~|~N2O4
N2O
C2H2
HC=-C-
-C=-C-
=#

"""
    EPPR78(components_or_groups;
    idealmodel = BasicIdeal,
    alpha = PR78Alpha,
    mixing = PPR78Rule,
    activity = nothing,
    translation = NoTranslation,
    userlocations = String[],
    ideal_userlocations = String[],
    alpha_userlocations = String[],
    mixing_userlocations = String[],
    translation_userlocations = String[],
    reference_state = nothing,
    verbose = false)
Enhanced Predictive Peng Robinson equation of state. it uses the following models:
- Translation Model: [`NoTranslation`](@ref)
- Alpha Model: [`PR78Alpha`](@ref)
- Mixing Rule Model: [`PPR78Rule`](@ref)
## References
1. Jaubert, J.-N., Privat, R., & Mutelet, F. (2010). Predicting the phase equilibria of synthetic petroleum fluids with the PPR78 approach. AIChE Journal. American Institute of Chemical Engineers, 56(12), 3225–3235. [doi:10.1002/aic.12232](https://doi.org/10.1002/aic.12232)
2. Jaubert, J.-N., Qian, J.-W., Lasala, S., & Privat, R. (2022). The impressive impact of including enthalpy and heat capacity of mixing data when parameterising equations of state. Application to the development of the E-PPR78 (Enhanced-Predictive-Peng-Robinson-78) model. Fluid Phase Equilibria, (113456), 113456. [doi:10.1016/j.fluid.2022.113456](https://doi.org/10.1016/j.fluid.2022.113456)
"""
function EPPR78(components;
    idealmodel = BasicIdeal,
    userlocations = String[],
    group_userlocations = String[],
    ideal_userlocations = String[],
    alpha_userlocations = String[],
    mixing_userlocations = String[],
    translation_userlocations = String[],
    reference_state = nothing,
    verbose = false)

    mixing = PPR78Rule(components,
            userlocations=mixing_userlocations,
            group_userlocations =group_userlocations,
            verbose = verbose)

    _components = mixing.groups.components
    alpha = PR78Alpha
    translation = NoTranslation

    return PR(_components;
    idealmodel = idealmodel,
    alpha = alpha,
    mixing = mixing,
    activity = nothing,
    translation = translation,
    userlocations = userlocations,
    ideal_userlocations = ideal_userlocations,
    alpha_userlocations = alpha_userlocations,
    mixing_userlocations = mixing_userlocations,
    translation_userlocations = translation_userlocations,
    reference_state = reference_state,
    verbose = verbose)
end
export EPPR78