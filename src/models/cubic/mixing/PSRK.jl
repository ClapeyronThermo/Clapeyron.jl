struct PSRKRule{γ} <: MHV1RuleModel
    components::Array{String,1}
    activity::γ
    references::Array{String,1}
end

@registermodel PSRKRule

"""
    PSRKRule{γ} <: MHV1RuleModel

    PSRKRule(components::Vector{String};
    activity = Wilson,
    userlocations::Vector{String}=String[],
    activity_userlocations::Vector{String}=String[],
    verbose::Bool=false)

## Input Parameters

None

## Input models

- `activity`: Activity Model

## Description

Mixing Rule used by the Predictive Soave-Redlich-Kwong [`PSRK`](@ref) EoS,
derived from the First Order modified Huron-Vidal Mixing Rule.
"""
PSRKRule

export PSRKRule
function PSRKRule(components::Vector{String}; activity = PSRKUNIFAC, userlocations::Vector{String}=String[],activity_userlocations::Vector{String}=String[], verbose::Bool=false)
    _activity = init_model(activity,components,activity_userlocations,verbose)
    references = String["10.1016/0378-3812(91)85038-V"]
    model = PSRKRule(components, _activity,references)
    return model
end

MHV1q(::PSRKRule,::RKModel) = 0.64663
