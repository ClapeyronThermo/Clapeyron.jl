# Extensions – CoolProp

`Clapeyron.jl` supports using all the single fluid helmholtz models available in [CoolProp](https://github.com/CoolProp/CoolProp), by loading the corresponding JSON data into a [`SingleFluid`](@ref) or [`MultiFluid`](@ref) model.

```julia
# let's try to load the R134a fluid model available in CoolProp
julia> model = SingleFluid("r134a")
ERROR: cannot found component file pentane. Try loading the CoolProp library (`using CoolProp`).

julia> using CoolProp  # the CoolProp extension is loaded automatically when CoolProp is loaded

julia> model = SingleFluid("r134a")
MultiParameter Equation of state for pentane:
 Polynomial power terms: 5
 Exponential terms: 6
 Gaussian bell-shaped terms: 5
```

The properties calculated by `Clapeyron.jl` match the ones calculated by `CoolProp`:

```julia
julia> mass_density(model,1e5,200.0)
1510.617461131427

julia> PropsSI("D","P",1e5,"T",200.0,"R134a")
1510.6174611314268
```

For ideal properties, the `SingleFluidIdeal` and `EmpiricIdeal` models are available.
Those models can be combined with another residual model:

```julia
# reference equations of state, ideal part
id = MultiFluidIdeal(["r134a","carbon dioxide"])
# PCSAFT for the mixture model
model = PCSAFT(["r134a","carbon dioxide"],idealmodel = id)
```
