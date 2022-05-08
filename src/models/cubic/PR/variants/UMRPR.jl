abstract type UMRPRModel <: PRModel end

UMRPR_SETUP = ModelOptions(
        :UMRPR;
        supertype=UMRPRModel,
        parent=PR_SETUP,
        members=[
            ModelMember(:alpha, :MTAlpha),
            ModelMember(:activity, :UMRRule),
            ModelMember(:mixing, :UNIFAC),
            ModelMember(:translation, :MTTranslation),
            ModelMember(:idealmodel, :BasicIdeal; groupcontribution_allowed=true),
        ],
        references=["10.1021/ie049580p"],
    )

createmodel(UMRPR_SETUP; verbose=true)
export UMRPR

"""
    UMRPR(components::Vector{String}; idealmodel=BasicIdeal,
    alpha = MTAlpha,
    mixing = UMRRule,
    activity = UNIFAC,
    translation=MTTranslation,
    userlocations=String[], 
    ideal_userlocations=String[],
    alpha_userlocations = String[],
    mixing_userlocations = String[],
    activity_userlocations = String[],
    translation_userlocations = String[],
    verbose=false)

Universal Mixing Rule Peng Robinson equation of state. it uses the following models:

- Translation Model: [`MTTranslation`](@ref)
- Alpha Model: [`MTAlpha`](@ref)
- Mixing Rule Model: [`UMRRule`](@ref) with (`UNIFAC`)(@ref) activity

## References

1. Voutsas, E., Magoulas, K., & Tassios, D. (2004). Universal mixing rule for cubic equations of state applicable to symmetric and asymmetric systems: Results with the Peng−Robinson equation of state. Industrial & Engineering Chemistry Research, 43(19), 6238–6246. doi:10.1021/ie049580p

"""
