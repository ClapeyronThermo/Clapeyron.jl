# v0.5.10

## New Features
- Association models don't allocate anymore in the case of a single association site pair.
- `saturation_pressure(model,T)` (`ChemPotVSaturation,IsoFugacitySaturation`) does not allocate if the calculation does not require a critical point calculation. Note that the function can still allocate if the EoS model itself allocates. the same optimizations were applied to `saturation_temperature` (`AntoineSaturation`,`ClapeyronSaturation`), `sublimation_pressure` and `melting_pressure`.
- Bulk properties now accept a `vol0` initial point for the volume solver.
- SAFT-VR-Mie uses a divided approach for calculating `d`: if θ = ℂ*ϵᵢ/T > 1, then it uses a 10-point gauss-laguerre integrator. Otherwise, the Assen method of finding a cut point and integrating the rest is used. A description of the method is found here: https://teqp.readthedocs.io/en/latest/models/SAFT-VR-Mie.html. the cut allows for better accuracy at higher reduced temperatures.

## Bug fixes
- Peng-Robinson now uses more accurate `Ωa` and `Ωb` values
- CPA/sCPA now uses SI units as input.
