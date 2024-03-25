# v0.5.11

## New Features
- Support for reference states. A reference state is a point in V-T space where H = H₀ and S = S₀. Setting those have uses in Reaction equilibria and when comparing between different models.In particular, Reference states are stored in a `ReferenceState <: ClapeyronParam` in the ideal model parameters. The `BasicIdeal` model is, intentionally, the only ideal model in Clapeyron that does not have this struct and, as a consequence, it is not able to set reference states.
- Support for superancillaries via [`EoSSuperancillaries.jl`](https://github.com/ClapeyronThermo/EoSSuperancillaries.jl). When the package is loaded, initial saturation points for cubics and PCSAFT are overloaded to use superancillary evaluations instead of the general `x0_sat_pure` function. in the case of `PCSAFT` models, it also speeds up the evaluation of `crit_pure`.
- New EoS: EOS-CG (2021) (`EOS_CG`), a reference model for humid gases and CCS mixtures.
- New EoS: Lee-Kesler-Plöcker (`LKP`)
- New EoS: Shomate ideal model (`ShomateIdeal`)
- database: PCPSAFT,gcPCSAFT and gcPCPSAFT are updated to use the values of Rehner (2023).
- new functions: `helmholtz_free_energy_res`,`gibbs_free_energy_res`,`internal_energy_res`, `enthalpy_res`
- database: `ReidIdeal` now uses the poling coefficients by default.
- database: `JobackIdeal` has support for more common group fragments used in gcPCSAFT.
- `melting_temperature`, `sublimation_temperature` does not allocate anymore. Note that the function can still allocate if the EoS model itself allocates.

## bug fixes
- Incorrect value for CPA with water (#256)
- Bug in SAFT-VR-SW (#165)
- Bug in CP-PC-SAFT
