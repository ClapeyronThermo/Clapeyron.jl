abstract type PR78Model <: PRModel end

PR78_SETUP = ModelOptions(
        :PR78;
        supertype=PR78Model,
        parent=PR_SETUP,
        members=[
            ModelMember(:alpha, :PR78Alpha),
            ModelMember(:activity, :Nothing),
            ModelMember(:mixing, :vdW1fRule),
            ModelMember(:translation, :NoTranslation),
            ModelMember(:idealmodel, :BasicIdeal; groupcontribution_allowed=true),
        ],
        references=["10.1021/I160057A011"],
    )

createmodel(PR78_SETUP; verbose=true)
export PR78

"""
    PR78(components::Vector{String}; idealmodel=BasicIdeal,
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

Peng Robinson (1978) equation of state. it uses the following models:

- Translation Model: [`NoTranslation`](@ref)
- Alpha Model: [`PR78Alpha`](@ref)
- Mixing Rule Model: [`vdW1fRule`](@ref)

## References

1. Robinson DB, Peng DY. The characterization of the heptanes and heavier fractions for the GPA Peng-Robinson programs. Tulsa: Gas Processors Association; 1978
"""
