struct PSRKRule{γ} <: MHV1RuleModel
    components::Array{String,1}
    activity::γ
    references::Array{String,1}
end

"""
    PSRKRule{γ} <: MHV1RuleModel

    PSRKRule(components;
    activity = Wilson,
    userlocations = String[],
    activity_userlocations = String[],
    verbose::Bool=false)

## Input Parameters

None

## Input models

- `activity`: Activity Model

## Description

Mixing Rule used by the Predictive Soave-Redlich-Kwong [`PSRK`](@ref) EoS,
derived from the First Order modified Huron-Vidal Mixing Rule.

## Model Construction Examples
```
# Note: this model was meant to be used exclusively with the PSRKUNIFAC activity model.

# Using the default database
mixing = PSRKRule(["water","carbon dioxide"]) #default: PSRKUNIFAC Activity Coefficient.
mixing = PSRKRule(["water","carbon dioxide"],activity = NRTL) #passing another Activity Coefficient Model
mixing = PSRKRule([("ethane",["CH3" => 2]),("butane",["CH2" => 2,"CH3" => 2])],activity = UNIFAC) #passing a GC Activity Coefficient Model.

# Passing a prebuilt model

act_model = NRTL(["water","ethanol"],userlocations = (a = [0.0 3.458; -0.801 0.0],b = [0.0 -586.1; 246.2 0.0], c = [0.0 0.3; 0.3 0.0]))
mixing = PSRKRule(["water","ethanol"],activity = act_model)

# Using user-provided parameters

# Passing files or folders
mixing = PSRKRule(["water","ethanol"]; activity = NRTL, activity_userlocations = ["path/to/my/db","nrtl_ge.csv"])

# Passing parameters directly
mixing = PSRKRule(["water","ethanol"];
                activity = NRTL,
                userlocations = (a = [0.0 3.458; -0.801 0.0],
                    b = [0.0 -586.1; 246.2 0.0],
                    c = [0.0 0.3; 0.3 0.0])
                )
```
"""
PSRKRule

export PSRKRule
function PSRKRule(components; activity = PSRKUNIFAC, userlocations = String[],activity_userlocations = String[], verbose::Bool=false)
    _activity = init_mixing_act(activity,components,activity_userlocations,verbose)
    references = String["10.1016/0378-3812(91)85038-V"]
    model = PSRKRule(format_components(components), _activity,references)
    return model
end

MHV1q(::PSRKRule,::RKModel) = 0.64663

function PSRKUNIFAC end
"""
    function PSRK(components;
    idealmodel = BasicIdeal,
    alpha = SoaveAlpha,
    mixing = PSRKRule,
    activity = PSRKUNIFAC,
    translation = PenelouxTranslation,
    userlocations = String[],
    ideal_userlocations = String[],
    alpha_userlocations = String[],
    mixing_userlocations = String[],
    activity_userlocations = String[],
    translation_userlocations = String[],
    reference_state = nothing,
    verbose = false)

## Description
Predictive Soave-Redlich-Kwong equation of state. it uses the following models:

- Translation Model: [`PenelouxTranslation`](@ref)
- Alpha Model: [`SoaveAlpha`](@ref)
- Mixing Rule Model: [`PSRKRule`](@ref) with [`PSRKUNIFAC`](@ref) activity model

##  References

1. Horstmann, S., Jabłoniec, A., Krafczyk, J., Fischer, K., & Gmehling, J. (2005). PSRK group contribution equation of state: comprehensive revision and extension IV, including critical constants and α-function parameters for 1000 components. Fluid Phase Equilibria, 227(2), 157–164. [doi:10.1016/j.fluid.2004.11.002](https://doi.org/10.1016/j.fluid.2004.11.002)
"""
function PSRK(components;
    idealmodel = BasicIdeal,
    userlocations = String[],
    group_userlocations = String[],
    ideal_userlocations = String[],
    alpha_userlocations = String[],
    mixing_userlocations = String[],
    activity_userlocations = String[],
    translation_userlocations = String[],
    reference_state = nothing,
    verbose = false)

    activity = PSRKUNIFAC(components,
    userlocations = activity_userlocations,
    group_userlocations = group_userlocations,
    verbose = verbose)

    _components = activity.groups.components #extract pure component list

    alpha = SoaveAlpha
    mixing = PSRKRule
    translation = PenelouxTranslation

    return RK(_components;
    idealmodel = idealmodel,
    alpha = alpha,
    mixing = mixing,
    activity = activity,
    translation = translation,
    userlocations = userlocations,
    ideal_userlocations = ideal_userlocations,
    alpha_userlocations = alpha_userlocations,
    mixing_userlocations = mixing_userlocations,
    activity_userlocations = activity_userlocations,
    translation_userlocations = translation_userlocations,
    reference_state = reference_state,
    verbose = verbose)
end
export PSRK
