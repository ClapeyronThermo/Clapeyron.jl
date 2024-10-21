```@meta
CurrentModule = Clapeyron
```

## Definitions

Within `Clapeyron`, we provide a few methods which we use to obtain thermodynamic properties (you can find out more details on how we implement these methods in our background information):

- `[volume](@ref)(model, p, T, z)`: Obtains the volume of a system at a given temperature, pressure and composition.
  If the phase is unknown, it will find the vapour and liquid roots and return the one that minimises the Gibbs free energy.
  This function is called by all of our bulk property methods.
- `[saturation_pressure](@ref)(model, T)`: Obtains the saturation pressures and volumes for a pure species.
- `[crit_pure](@ref)(model)`: Obtains the critical point for a pure species.

Clearly this is not an exhaustive list and we make absolutely no guarantees as to the quality of these methods.
However, if you have a new method or algorithm you would like to implement or test out, it is possible to do so.

## Custom initial guesses example

One of the most common reasons for our methods to fail (that we've seen) is due to poor initial guesses.
If you are experiencing issues with our `saturation_pressure` method specifically, you could try modifying the initial guess for a particular equation of state (see the custom models documentation for explanation on abstract types):

```julia
function Clapeyron.x0_sat_pure(model::PCSAFTModel,T,z=SA[1.0])
  # Obtain the volume lower bound for that particular system
  Vlb = lb_volume(model,z)*one(T)

  # Relative to the lower bound, define your initial guesses.
  # We log10 the results as our solvers solve for the log10 of the volume.
  return Vlb*1.5,Vlb*100
end
```

Note that we do need prefix [`x0_sat_pure`](@ref) with `Clapeyron.` as we do not export this function normally; including this function in our script with force Clapeyron to use it instead of the default.
These modifications can also be done for `x0_crit_pure` although `crit_pure` has proven to be quite reliable.

## Custom volume solver

For something a bit more substantial, you can also modify the volume function itself.
We do this with cubics (there is a direct procedure for obtaining the volume roots) and ideal models (the volume is the ideal gas).
The default solver method ([`volume`](@ref)) can be overloaded by defining:

```julia
function Clapeyron.volume_impl(model::MyEoSModel,p,T,z=SA[1.0],phase=:unknown,threaded=false,vol0=nothing)

  # INSERT YOUR ALGORITHM HERE

  return vol
end
```

Clapeyron will automatically call your implementation when your model is evaluated.

## Custom saturation solver

For saturation solvers ([`saturation_pressure`](@ref),[`saturation_temperature`](@ref)) you can dispatch on a different saturation method.
Let's create one, that just evaluates Antoine coefficients:

```julia
struct DirectAntoine{C} <: Clapeyron.SaturationMethod
  A::Float64
  B::Float64
  C::Float64
  crit::C
end

# Defining interface, check Clapeyron.SaturationMethod for more information.
# All saturation methods require passing crit as a keyword.

DirectAntoine(A,B,C;crit = nothing) = DirectAntoine(A,B,C,crit)

function Clapeyron.saturation_temperature_impl(model::EoSModel,T,method::DirectAntoine)
  A = method.A
  B = method.B
  C = method.C
  p = exp(A + B/(T+C))
  vl = volume(model,p,T,phase = :l)
  vv = volume(model,p,T,phase = :v)
  return p,vl,vv
end
```

You can now call `saturation_pressure(model,T,method)` where `method = DirectAntoine(A,B,C)`.
At the moment, the default algorithm iterates directly over volumes ([`ChemPotVSaturation`](@ref)), but we also have saturation via isofugacity ([`IsofugacitySaturation`](@ref)) and superancillaries for cubics ([`SuperAncSaturation`](@ref)).
The same thing can be done with `saturation_temperature`.

## Custom TP-Flash solver

We support the same procedure to define your custom Temperature-pressure flash methods:

```julia
struct MyRachfordRice{K} <: Clapeyron.TPFlashMethod
  K0::K
end

# defining interface, check Clapeyron.TPFlashMethod for more information
numphases(::MyRachfordRice) = 2
# we perform index reduction, to create smaller models in case one component has zero composition.
function Clapeyron.index_reduction(method::MyRachfordRice,non_zero_indices)
  return MyRachfordRice(index_reduction(method.K0,non_zero_indices))
end

function Clapeyron.tp_flash_impl(model::EoSModel,p,T,z,method::MyRachfordRice)
  # perform Rachford-Rice,returns x, y, α₀
  # ...
  X = vcat(x',y')
  n = X.*[1-α₀
                α₀]  .* sum(z)
  g = (gibbs_free_energy(model,p,T,x)*(1-α₀)+gibbs_free_energy(model,p,T,y)*α₀)/R̄/T
  return X,n,g
end
```

## I have a better method...

If your custom methods end up being more-efficient than ours or you develop one that we do not currently support, please do start a pull request and we will gladly add it to the package!
