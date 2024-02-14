"""
    HelmAct::MultiFluid
    HelmAct(components;
        pure_userlocations = String[],
        activity = PSRKUNIFAC,
        activity_userlocations = String[],
        estimate_pure = false,
        coolprop_userlocations = false,
        Rgas = Rgas,
        reference_state = nothing,
        verbose = verbose)

## input Parameters

- CoolProp JSON pure fluid files

## Input Models
- `activity`: activity model.

## Description

Creates a Helmholtz + Activity (HelmAct) model:

```
aᵣ = ∑xᵢaᵣᵢ(δ,τ) + Δa
Δa = gᴱᵣ/RT - log(1+bρ)/log(1+bρref) * ∑xᵢ(aᵣᵢ(δref,τ) - aᵣᵢ(δrefᵢ,τᵢ))
τᵢ = Tcᵢ/T
δref = ρref/ρr
δrefᵢ = ρrefᵢ/ρcᵢ
b = 1/1.17ρref
```

where `gᴱᵣ` is the residual part of the excess gibbs free energy obtained from an activity model.

## References
1. Jäger, A., Breitkopf, C., & Richter, M. (2021). The representation of cross second virial coefficients by multifluid mixture models and other equations of state. Industrial & Engineering Chemistry Research, 60(25), 9286–9295. [doi:10.1021/acs.iecr.1c01186](https://doi.org/10.1021/acs.iecr.1c01186)


"""
function HelmAct(components;
    pure_userlocations = String[],
    activity = PSRKUNIFAC,
    activity_userlocations = String[],
    idealmodel = nothing,
    ideal_userlocations = String[],
    estimate_pure = false,
    coolprop_userlocations = false,
    Rgas = Rgas,
    reference_state = nothing,
    verbose = verbose)

    init_activity = init_model(activity,components,activity_userlocations,verbose)
    if has_groups(init_activity)
        components = init_activity.groups.components
    else
        components = init_activity.components
    end
    departure = GEDeparture(init_activity)
    return MultiFluid(components;
        pure_userlocations = pure_userlocations,
        departure = departure,
        mixing = LinearMixing(),
        estimate_pure = estimate_pure,
        coolprop_userlocations = coolprop_userlocations,
        Rgas = Rgas,
        reference_state = reference_state,
        verbose = false,
        idealmodel = idealmodel,
        ideal_userlocations = ideal_userlocations,
        )
end

export HelmAct