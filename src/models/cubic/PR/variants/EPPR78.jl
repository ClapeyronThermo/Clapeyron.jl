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
    EPPR78(components::Vector{String}; idealmodel=BasicIdeal,
    alpha = PR78Alpha,
    mixing = vdW1fRule,
    activity = nothing,
    translation=NoTranslation,
    userlocations=String[], 
    ideal_userlocations=String[],
    alpha_userlocations = String[],
    mixing_userlocations = String[],
    translation_userlocations = String[],
    verbose=false)

Enhanced Predictive Peng Robinson equation of state. it uses the following models:

- Translation Model: `NoTranslation`
- Alpha Model: `PR78Alpha`
- Mixing Rule Model: `PPR78Rule`

## References

1. Robinson DB, Peng DY. The characterization of the heptanes and heavier fractions for the GPA Peng-Robinson programs. Tulsa: Gas Processors Association; 1978
"""
function EPPR78(components::Vector{String}; idealmodel=BasicIdeal,
    alpha = PR78Alpha,
    mixing = PPR78Rule,
    activity = nothing,
    translation=NoTranslation,
    userlocations=String[], 
    ideal_userlocations=String[],
    alpha_userlocations = String[],
    mixing_userlocations = String[],
    translation_userlocations = String[],
    verbose=false)

    return PR(components;
    idealmodel = idealmodel,
    alpha = alpha,
    mixing=mixing,
    activity = activity,
    translation=translation,
    userlocations = userlocations,
    ideal_userlocations = ideal_userlocations,
    alpha_userlocations = alpha_userlocations,
    mixing_userlocations = mixing_userlocations,
    translation_userlocations = translation_userlocations,
    verbose = verbose)
end
export EPPR78