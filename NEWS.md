# v0.3.12

## new features
- new EoS: RKPR (`RKPR`)
- new functions: `cross_second_virial` and `equivol_cross_second_virial`, to calculate B12, at specified z conditions, or by setting equivolumetric mixing.
- cubic EoS now support `recombine!`
- `PenelouxTranslation` and `RackettTranslation` now cache their results, resulting in a speed up when using volume translated EoS.
- `Clapeyron.cite` now accepts the optional argument `out`, that can be `doi` (current default, DOI) or `:bib` (for BibTeX)
- 
## Bug fixes
- various bug fixes in SAFT-VRQ-Mie. the Hard sphere term now performs the non-additive mixing rule for the multicomponent case. initializations with integer values are fixed.
