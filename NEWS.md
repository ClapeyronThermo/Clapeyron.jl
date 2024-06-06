# v0.5.12

## New Features
- New method: `MultiPhaseTPFlash`, that solves multiphase,multicomponent TP-flash with automatic phase detection. this method is now the default when calling `tp_flash`
- New method: `Tproperty(model,p,prob,z,property)` to calculate temperatures, given pressure and another property.
- additional method: `x0_volume_liquid(model,p,T,z)` and `x0_volume_solid(model,p,T,z)` can be overloaded to calculate liquid an solid volumes, using the pressure as information. They are defined as `x0_volume_liquid(model,p,T,z) = x0_volume_liquid(model,T,z)` and `x0_volume_solid(model,p,T,z) = x0_volume_solid(model,T,z)`
- tangent plane distance (`tpd`) calculations are now faster.
- `VT_diffusive_stability` now uses `eigmin` instead of the full eigen calculation.
- `isstable` now works on (P,T,z) space, for the (V,T,z) space, use `VT_isstable`. there are now (P,T,z) versions of each stability function.
- calculation of volumes,saturation pressures and critical points of CPA models now defaults to the inner cubic model when there is no association present.
- The default association implementation now uses a combination of accelerated successive substitution and newton optimization. While increasing allocations, the method is faster.
- the default `volume` implementation now uses implicit AD to support derivatives. instead of propagating derivative information through the iterative procedure. This allows workloads of the type: `ForwardDiff.derivative(_p -> property(model,_p,T,z,phase = :l,vol0 = v0),p)` to be efficiently calculated.

## Bug fixes
- PCPSAFT: typo in unlike asssociation parameters
