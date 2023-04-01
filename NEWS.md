# v0.4.9

## New Features

- `ideal_consistency(model,V,T,z)` that checks if `da0/dV + sum(z)/V` is zero (or as close to zero as the Floating Point Format allows it.)
- broadcasting on `AssocParam` is defined (`bondvol .= 1e-10 .* bondvol .^3`) 

## Bug Fixes

- proper namespaces in `@registermodel` ([#161](https://github.com/ClapeyronThermo/Clapeyron.jl/issues/161))
- fixed bug in `MichelsenTPFlash` when using non-volatiles and `second_order = false`
- fixed bug when building `UNIFACFVPoly`