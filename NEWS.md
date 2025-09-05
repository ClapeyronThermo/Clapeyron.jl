# v0.6.16

## New Features

- New model type: Gibbs models (`model <: GibbsBasedModel`): EoS that use the gibbs energy as the main function. Properties and single component equilibria are implemented.
- Liquid Volume models are now a subtype of `GibbsBasedModel`
- New method: `component_list(model)` that returns a vector of component names.
- New property method: `sle_solubility_T`
- New Gibbs EoS: IAPWS-06 standard for ice-ih (`IAPWS06`)
- New Gibbs EoS: Grenke-Elliott model for liquid water (`GrenkeElliottWater`)
- New Gibbs EoS: Holten 2-state model for liquid supercooled water (`HoltenWater`)
- New Gibbs EoS: Jäger and Span model for solid carbon dioxide (`JagerSpanSolidCO2`)
- New Activity model: Van Laar binary activity coefficient model (`VanLaar`)
- New Activity model: Margules binary activity coefficient model (`Margules`)
- New Activity model: Flory-Huggins activity coefficient model (`FloryHuggins`)
- New Correlation Model: Antoine saturation pressure correlation (`AntoineEqSat`)
- Reference State: new type: `:ideal_gas`, that calculated the reference state only with the ideal gas model.
- Misc documentation improvements.

## Bug Fixes

- Volume solver: the solver now tries to detect if there was a "jump" in the chosen (liquid or vapour) branch being iterated. For EoS with more than two maxwell loops, branch jumps sometimes caused the final volume to be in a invalid branch.
- Misc improvements to the database.
- Fixed initialization for `ReferenceState` when `H0` and `S0` are vectors.
- improved initialization for electrolyte models.
- bubble/dew: improved handling of negative fractions.
