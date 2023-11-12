# v0.5.8

## New Features
- `Base.getindex` and `Base.setindex` with `SingleParam`, `PairParam` and `AssocParam` now works with strings. the strings are compared with the components (or groups) stored in each param. in particular `AssocParam` allows set/get index methods if you pass a `Tuple{String,String}`:
```julia
julia> model = PPCSAFT(["water","ethanol"],assoc_options = AssocOptions(combining = :esd))
PPCSAFT{BasicIdeal} with 2 components:
 "water"
 "ethanol"
Contains parameters: Mw, segment, sigma, epsilon, dipole, dipole2, epsilon_assoc, bondvol

julia> model.params.bondvol[("water","a"),("water","b")]
0.35319

julia> model.params.bondvol[("water","a"),("water","b")] = 0.36
0.36

julia> model.params.bondvol[("water","a"),("water","b")]
0.36
```
- `PCPSAFT` is defined (alias for `PPCSAFT`)
- New EOS: Critical-point based PC-SAFT `CPPCSAFT` (https://doi.org/10.1021/ie502633e)

## Bug Fixes
- bug in ether and aldehyde parameters in UNIFAC (https://github.com/ClapeyronThermo/Clapeyron.jl/issues/225)

