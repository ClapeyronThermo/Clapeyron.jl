# v0.4.13

## New Features

- (Experimental) Initial `Symbolics.jl` support for bulk properties on `EoSModel`. In particular, all Activity models, cubic models, SAFT models without association, and empiric models are supported. for example, this is now supported:

```julia
@variables y(..)[1:4], n(..)[1:4]
mixture = UNIFAC(["water","ethanol","methanol","1-propanol"])
moles = [n(t)[i] for i in 1:4] #Clapeyron accepts mole amounts, so it is not necessary to perform transformations to mole fractions
bc_l = y(t, 0.0) .~ activity_coefficient(mixture, 1.0, 298.15, moles)
```

## Bug Fixes

- automatic precompile is disabled in this version.