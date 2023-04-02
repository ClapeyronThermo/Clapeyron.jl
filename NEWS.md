# v0.4.9

## New Features

- `ideal_consistency(model,V,T,z)` that checks if `da0/dV + sum(z)/V` is zero (or as close to zero as the Floating Point Format allows it.)
- broadcasting on `AssocParam` is defined (`bondvol .= 1e-10 .* bondvol .^3`) 
- you can pass functions that build models instead of  EoSModel types. for example, now you can do:
    ```julia
    function myvdW(components;userlocations = String[],verbose = false)
        return vdW(components;userlocations = userlocations,verbose = verbose,alpha = SoaveAlpha)
    end

    model = Wilson(["water","ethanol"];puremodel=myvdW)
    ```

## Bug Fixes

- proper namespaces in `@registermodel` ([#161](https://github.com/ClapeyronThermo/Clapeyron.jl/issues/161))
- fixed bug in `MichelsenTPFlash` when using non-volatiles and `second_order = false`
- fixed bug when building `UNIFACFVPoly`
