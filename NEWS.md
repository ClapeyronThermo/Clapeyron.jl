# v0.6.0

## New Features
- New models: Electrolyte models are now supported! We have introduced the `ESElectrolyte` framework which will let users combine any electrostatic model (`DH`, `MSA` and `Born`) and relative static permittivity model with any of our supported equations of state. Due to this flexibility, we now support four existing SAFT-type electrolyte equations (with planned support for more):
  - `ePCSAFT`
  - `SAFTVREMie`
  - `eSAFTVRMie`
  - `SAFTgammaEMie`
- New method: Two new methods specific to electrolytes have been added: `mean_ionic_activity_coefficient` and `osmotic_coefficient`, along with their saturated variants.
- New method: `MultiPhaseTPFlash`, that solves multiphase,multicomponent TP-flash with automatic phase detection. this method is now the default when calling `tp_flash` with more than two components and helmholtz-based models.
- New method: `Tproperty(model,p,prob,z,property)` to calculate temperatures, given pressure and another property.
- New model: To model solubility of salts, `SolidKs` has been added in order to obtain the solubility of salts using the infinite-dilution approach as opposed to the pure-fluid approach using `SolidHfus`.
- additional method: `x0_volume_liquid(model,p,T,z)` and `x0_volume_solid(model,p,T,z)` can be overloaded to calculate liquid an solid volumes, using the pressure as information. They are defined as `x0_volume_liquid(model,p,T,z) = x0_volume_liquid(model,T,z)` and `x0_volume_solid(model,p,T,z) = x0_volume_solid(model,T,z)`
- tangent plane distance (`tpd`) calculations are now faster.
- `VT_diffusive_stability` now uses `eigmin` instead of the full eigen calculation.
- `isstable` now works on (P,T,z) space, for the (V,T,z) space, use `VT_isstable`. there are now (P,T,z) versions of each stability function.
- calculation of volumes,saturation pressures and critical points of CPA models now defaults to the inner cubic model when there is no association present.
- The default association implementation now uses a combination of accelerated successive substitution and newton optimization. While increasing allocations, the method is faster.
- the default `volume` implementation now uses implicit AD to support derivatives. instead of propagating derivative information through the iterative procedure. This allows workloads of the type: `ForwardDiff.derivative(_p -> property(model,_p,T,z,phase = :l,vol0 = v0),p)` to be efficiently calculated.
- `Clapeyron.tpd` code has been optimized. `tpd` has new keywords: `break_first`, that tries to return a negative tpd as early as possible, `lle` for only calculating TPD in liquid phases, `strategy`, that changes the search strategy between a K-value search (`:wilson`), a pure component search (`:pure`) or both strategies (`default`).
- `Clapeyron.tpd` now supports activity models (if the keyword `lle` is set to `true`)
- New EoS: modified Lee-Kesler-Plöcker with consistent parameters (`LKPmod`)
- New EoS: Lee-Kesler-Plöker-equation of state, Sabozin-Jäger-Thol enhancement (`LKPSJT`, `enhancedLKP`)
## Bug fixes
- PCPSAFT: typo in unlike asssociation parameters
