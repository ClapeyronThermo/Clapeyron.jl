# v0.5.12

## New Features
- New method: `Tproperty(model,p,prob,z,property)` to calculate temperatures, given pressure and another property.
- additional method: `x0_volume_liquid(model,p,T,z)` and `x0_volume_solid(model,p,T,z)` can be overloaded to calculate liquid an solid volumes, using the pressure as information. They are defined as `x0_volume_liquid(model,p,T,z) = x0_volume_liquid(model,T,z)` and `x0_volume_solid(model,p,T,z) = x0_volume_solid(model,T,z)`