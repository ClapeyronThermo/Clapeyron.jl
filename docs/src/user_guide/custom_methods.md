## Definitions

Within `Clapeyron`, we provide a few methods which we use to obtain thermodynamic properties (you can find out more details on how we implement these methods in our background information):

- `volume(model, p, T, z)`: Obtains the volume of a system at a given temperature, pressure and composition. If the phase is unknown, it will find the vapour and liquid roots and return the one that minimises the Gibbs free energy. This function is called by all of our bulk property methods.
- `saturation_pressure(model, T)`: Obtains the saturation pressures and volumes for a pure species.
- `crit_pure(model)`: Obtains the critical point for a pure species.

Clearly this is not an exhaustive list and we make absolutely no guarantees as to the quality of these methods. However, if you have a new method or algorithm you would like to implement or test out, it is possible to do so.

## Custom initial guesses example

One of the most common reasons for our methods to fail (that we've seen) is due to poor initial guesses. If you are experiencing issues with our `sat_pure` method specifically, you could try modifying the initial guess for a particular equation of state (see the custom models documentation for explanation on abstract types):

```julia
function Clapeyron.x0_sat_pure(model::PCSAFTModel,T,z=SA[1.0])
  # Obtain the volume lower bound for that particular system
  Vlb = lb_volume(model,z)*one(T)
  
  # Relative to the lower bound, define your initial guesses. We log10 the results as our solvers solve for the log10 of the volume.
  return log10.([Vlb*1.5,Vlb*100])
end
```

Note that we do need prefix `x0_sat_pure` with `Clapeyron.` as we do not export this function normally; including this function in our script with force Clapeyron to use it instead of the default. These modifications can also be done for `x0_crit_pure` although `crit_pure` has proven to be quite reliable.

## Custom volume solver example

For something a bit more substantial, you can also modify the volume function itself. Since this function is exported in Clapeyron, you do not need to prefix with `Clapeyron.`:

```julia
function volume(model::EoSModel,p,T,z=SA[1.0];phase=:unknown,threaded=true)
  
  # INSERT YOUR ALGORITHM HERE
  
  return vol
end
```

Clapeyron will automatically overwrite the default function and use this one instead. 

If your custom methods end up being more-efficient than ours or you develop one that we do not currently support, please do start a pull request and we will gladly add it to the package!
