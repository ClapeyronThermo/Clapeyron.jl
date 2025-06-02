# v0.6.13

## New Features

- Combining rules: support for matrices for all inplace combining rules
- improved speed for in bubbledew calculations with nonvolatiles/noncondensables when one phase has only one element. 
- new method: calculation of internal energy - volume flash (`uv_flash`) for single component models.
- `MultiFluid`: support for Double-Exponential terms.
- Cubics: new function, `CubicModel(modeltype,params,components;kwargs...)` that constructs a cubic model. All Clapeyron cubic models,as well as CPA, now use this function for their constructors.
- `SAFTgammaMie` new method: `SAFTgammaMie(groups::GroupParam,param::Dict{String,ClapeyronParam})` for easier construction.

## Bug Fixes

- CPA: improved recombining rules
- `SingleFluid`: derivatives of non-analytical terms at exactly the critical point are now calculated at an inifinitesimal point `(τ + εₜ,δ + εᵥ)`. This returns in somewhat better results than setting the whole term to zero.
- improved volume initial points for `SingleFluid`
- fixed bugs in `recombine!(model::SAFTgammaMie)`
- fixed bugs in `recombine!(model::SAFTVRMie)`
- Improvement of XY-flash results when one of the bubbledew calculation fails.
- fix bug when using MultiParameter EoS with tp-flash.
- assoc views: the sizes of assoc view now respect symmetry (diagonal entries have a square size) and transpose (`assoc_param[i,j] == transpose(assoc_param[j,i])`). Index access is not affected.
- 
