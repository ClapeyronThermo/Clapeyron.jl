## Definitions

Within `Clapeyron`, it is possible to create your own model/equation of state whilst still using all of the property estimation tools we provide. In order to do this, you must create the model. We provide a macro to make it easy to set it up. If you wish to create a new model called CustomEoS, you just need call the `@newmodel` macro with three parameters:

```julia
@newmodel CustomEoS CustomEosModel CustomEosParam
```

1. Model struct name
2. Model abstract type
3. Model parameters struct

We will now give a brief overview of each of these aspects.

### Model struct name

This is the concrete type, which is simply the name of the equation of state which will be used to generate the model. Within `Clapeyron`, we've tried to keep these names as obvious as possible (_e.g._ `vdW`, `PR`, `PCSAFT`, `SAFTVRMie`, `GERG2008`, _etc._). For your own model, this can be whatever you want as long is the identification makes sense to you.

### Model abstract type

In Julia, objects can only be sub-types of abstract types (and not concrete types), which will allow you to inherit the behaviour of the parent(s). In order to maintain a hierarchy of models, we highly encourage you to create an abstract type for your model that is a subtype one of these families of EoS, and dispatch on this newly created abstract type instead of on your model's concrete type. This is not enforced, but we suffix all abstract types in this category with `Model` as a convention. The abstract type that is at the top of the hierarchy is `EoSModel`. From this parent, we branch into more-specific EoS sub-types:

- `SAFTModel`: These are the models which have three parameters in common: segment size, $\sigma$, potential depth, $\epsilon$, and number of segments, $m$. All other SAFT-type models branch from this parent (`PCSAFTModel`, `SAFTVRMieModel`, `softSAFTModel`, _etc._). 
- `CubicModel`: These are the models whose parameters can be obtained from the critical temperature and pressure. With the exception of `CPAModel`, all cubics have a common structure where one can re-arrange the equation for the pressure as a third-order polynomial. As such, we define a subtype of `CubicModel`, `ABCubicModel` (_e.g._ `vdWModel`, `RKModel`, `SRKModel`, `PRModel`).
- `EmpiricHelmholtzModel`: These are the high-accuracy, multi-parameter models for specific species or systems (_e.g._ `GERG2008Model`, `IAPWS95Model`). There is no general structure to the models and they are treated as self-contained.
- `IdealModel`: Often overlooked, these models supplement the `SAFTModel` and `CubicModel` by providing the ideal contribution. Whilst the parameters and structure aren't usually the same between ideal models, this is unnecessary as the equation for the pressure is always $pV=Nk_\mathrm{B}T$ . 

For example, if we wish to create a new EoS model called `CustomEoS`, we will create an abstract type `CustomEoSModel`, that is a sub-type of another abstract type, say `CubicModel` (but it could also inherit from something further down the hierarchy, like `vdWModel`):

```
abstract type CustomEoSModel <: CubicModel end
```

For models that are sub-types of `SAFTModel` or `CubicModel`, most methods will be instantly compatible because methods will be able to make use of a set of the EoS parameters for finding the initial guesses. If your model is not of either of these types, it will be necessary to define a few additional functions:

- `lb_volume(model::CustomEoSModel,T,z)`: This must output the smallest possible value of the volume for your particular model using the model parameters. In SAFT equations, this is equivalent to a packing fraction of one ($\eta=1$) whilst in cubics, it is equivalent to the `b` parameter.
- `T_scale(model::CustomEoSModel,z)`: This must output the temperature scaling for your model using the model parameters. In SAFT equations, this is usually the potential depth whilst in cubics it is the critical temperature.
- `p_scale(model::CustomEoSModel,z)`: This must output the pressure scaling for your model using the model parameters. In cubics, this is the critical pressure whilst in SAFT we use a more complicated definition using the segment sizes and potential depths.

### Model parameters struct

The parameters for a particular system are all stored within a struct that is a subtype of `EoSParam`. By convention, we suffix these with `Param`. These structs should contain the model parameters, which comprise objects of types `SingleParam{T}`, `PairParam{T}`, and `AssocParam{T}`, where `T` is usually a base type (`Float64`, `Integer`, `String`, etc). Below is an example of a generic param struct for a SAFT and cubic model.

```julia
struct GenericSAFTParam <: EoSParam
    Mw::SingleParam{Float64}
    segment::SingleParam{Float64}
    sigma::PairParam{Float64}
    epsilon::PairParam{Float64}
    epsilon_assoc::AssocParam{Float64}
    bondvol::AssocParam{Float64}
end
```

```julia
struct GenericCubicParam <: EoSParam
    Tc::SingleParam{Float64}
    pc::SingleParam{Float64}
    Mw::SingleParam{Float64}
    a::PairParam{Float64}
    b::PairParam{Float64}
end
```

The parameters are wrapped in one of the following structs:

- `SingleParam{T}`: These are the parameters associated with a pure species or indexed by a single index ($i$). For example, the molecular weight, number of segments, critical temperature and pressure.
- `PairParam{T}`: These are the parameters associated with a pair of species or indexed by two indices ($i$ and $j$). If the two indices are the same ($i=j$), they refer to the like species, otherwise ($i\neq j$), they refer to the unlike species. For example, the cubic `a` and `b` parameters, the SAFT segment size (`sigma`) and potential depth (`epsilon`) parameters.
- `AssocParam{T}`: There are the parameters associated with both a pair of species and association sites (see background documentation for what these are). They are indexed by four indices (the species $i$ and $j$, and the sites $a$ and $b$). For example, potential depth of the association interaction (`epsilon_assoc`) and the bonding volume (`bondvol`).

If there exists a model parameter struct that contains exactly the same parameters as the one that you are about to create, you can also directly use that existing struct.

## PC-SAFT Example

Once all the above has been defined, we are ready to build our own model. Let us imagine we are trying to implement `PCSAFT`. 

1. We first define the Model name, type and parameters:

```julia
# Defining an abstract type for this model type
abstract type PCSAFTModel <: SAFTModel end

# Defining the parameters used by the model
struct PCSAFTParam <: EoSParam
    Mw::SingleParam{Float64}
    segment::SingleParam{Float64}
    sigma::PairParam{Float64}
    epsilon::PairParam{Float64}
    epsilon_assoc::AssocParam{Float64}
    bondvol::AssocParam{Float64}
end

# Creating a model struct called PCSAFT, which is a sub-type of PCSAFTModel, and uses parameters defined in PCSAFTParam
@newmodel PCSAFT PCSAFTModel PCSAFTParam
```

2. The next step is to create an outer constructor for the model that you have just defined. It should have the same name as the struct above. It can take the following arguments, although these can be hard-coded if you so wish:
   * `components`: A list of strings that identify the components.
   * `idealmodel`: An ideal model, but this can be left as an optional parameter.
   * `userlocations`: A list of strings that are paths to the databases that you are using.
   * `ideal_userlocations`: Same as above, but for ideal models.
   * `verbose`: For when you want to print more information to the console.
   * `assoc_options`: For use in the association sites solver.

```Julia
function PCSAFT(components; idealmodel=BasicIdeal, userlocations=String[], ideal_userlocations=String[], verbose=false,assoc_options = AssocOptions())
  	# Obtain a Dict of parameters. We pass in custom locations through the optional parameter userlocations.
    params,sites = getparams(components; userlocations=userlocations, verbose=verbose)
  
    # For clarity, we assign the contents of the returned dict to their own variables.
    segment = params["segment"]
    k = get(params,"k",nothing) #if k is not provided, it will be not be considered
    Mw = params["Mw"]
    # Here, we modify the values of the sigma parameter first.
    params["sigma"].values .*= 1E-10
  
    # In some cases, we may not have the unlike parameters and will need to use combining rules. You can also define your own combining rules for this.
    sigma = sigma_LorentzBerthelot(params["sigma"])
    epsilon = epsilon_LorentzBerthelot(params["epsilon"], k)
  
    epsilon_assoc = params["epsilon_assoc"]
    bondvol = params["bondvol"]
  
    bondvol,epsilon_assoc = assoc_mix(bondvol,epsilon_assoc,sigma,assoc_options) #combining rules for association. if you want to perform cross-association mixing, check the AssocOptions docs

    # Now we can create the parameter struct that we have defined.
    packagedparams = PCSAFTParam(Mw, segment, sigma, epsilon, epsilon_assoc, bondvol)
  
    # Although optional, it's generally good practise to cite your models!
    references = ["10.1021/ie0003887", "10.1021/ie010954d"]

    # Build the model.
    model = PCSAFT(packagedparams, sites, idealmodel; ideal_userlocations=ideal_userlocations, references=references, verbose=verbose,assoc_options = assoc_options)
  
    # Return the PCSAFT object that you have just created.
    return model
end
```

3. Define all of the model equations. We encourage you to use the full range of Unicode characters where it makes your code clearer to read!

   As convention, the first four arguments should be `model`, `V`, `T` and `z`; any other variables should come after.

   If we obey that convention, we may use the `@f` macro, which automatically substitutes the first four parameters for compactness. For example, `@f(func,i,j)` is equivalent to calling `func(model,V,T,z,i,j)`.

    Clapeyron obtains all the properties of a model by differenciating the total helmoltz energy ([`eos`](@ref)) or the residual helmoltz energy ([`eos_res`](@ref)).  `eos` and `eos_res` themselves are defined in terms of the reduced ideal helmholtz energy ([`a_res`](@ref)). In this case, we are going to define `a_res` for our own model:

   ```julia
   function Clapeyron.a_res(model::PCSAFTModel, V, T, z)
       return @f(a_hc) + @f(a_disp) + @f(a_assoc)
   end
   
   function a_hc(model::PCSAFTModel, V, T, z)
       x = z/∑(z)
       m = model.params.segment.values
       m̄ = ∑(x .* m)
       return m̄*@f(a_hs) - ∑(x[i]*(m[i]-1)*log(@f(g_hs,i,i)) for i ∈ @comps)
   end
   
   function d(model::PCSAFTModel, V, T, z, i)
       ϵii = model.params.epsilon.values[i,i]
       σii = model.params.sigma.values[i,i]
       return σii * (1 - 0.12exp(-3ϵii/T))
   end
   
   function ζ(model::PCSAFTModel, V, T, z, n)
       ∑z = ∑(z)
       x = z * (one(∑z)/∑z)
       m = model.params.segment.values
       res = N_A*∑z*π/6/V * ∑((x[i]*m[i]*@f(d,i)^n for i ∈ @comps))
   end
   
   function g_hs(model::PCSAFTModel, V, T, z, i, j)
       di = @f(d,i)
       dj = @f(d,j)
       ζ2 = @f(ζ,2)
       ζ3 = @f(ζ,3)
       return 1/(1-ζ3) + di*dj/(di+dj)*3ζ2/(1-ζ3)^2 + (di*dj/(di+dj))^2*2ζ2^2/(1-ζ3)^3
   end
   
   function a_hs(model::PCSAFTModel, V, T, z)
       ζ0 = @f(ζ,0)
       ζ1 = @f(ζ,1)
       ζ2 = @f(ζ,2)
       ζ3 = @f(ζ,3)
       return 1/ζ0 * (3ζ1*ζ2/(1-ζ3) + ζ2^3/(ζ3*(1-ζ3)^2) + (ζ2^3/ζ3^2-ζ0)*log(1-ζ3))
   end
   
   # INSERT REST OF CODE
   ```

4. With all the above defined in a single script, we can save the file as `PCSAFT.jl` and then include it in our jupyter notebooks (for example) and use the model with all of our existing method:

   ```julia
   include("PCSAFT.jl")
   
   model = PCSAFT(["carbon dioxide"])
   
   p = 20e6
   T = range(290,460,length=200)
   
   Cp = isobaric_heat_capacity.(model, p, T)
   
   (T_c, p_c, V_c) = crit_pure(model)
   
   T_sat = range(220,T_c,length=200)
   
   (p_sat, V_l_sat, V_v_sat) = saturation_pressure(model,T_sat)
   ```

## sPC-SAFT Example

Instead of developing an entirely new model, some of us may want to modify or extend an existing one. `sPCSAFT` is an example where we want to modify parts of regular `PCSAFT` but keep the rest the same. We can do this in a very succinct way making this new model a sub-type of the abstract type associated with another model.

- When we define the model type, `sPCSAFT` is a sub-type of `PCSAFT`:

  ```julia
  abstract type sPCSAFTModel <: PCSAFTModel end
  ```

- Since the parameters are the same, we can just use the same model params when creating the model:

  ```julia
  @newmodel sPCSAFT sPCSAFTModel PCSAFTParam
  ```

  This may not be the case if we're extending a model (_e.g._ if we're adding polar or ionic terms, we may need to define a new parameter struct to include the new parameters).

- When defining the model equations, we only need to write those that have been changed in `sPCSAFT`:

  ```julia
  function a_hc(model::sPCSAFTModel, V, T, z)
      x = z/sum(z)
      m = model.params.segment.values
      m̄ = ∑(x .* m)
      return m̄*@f(a_hs) - (m̄-1)*log(@f(g_hs))
  end
  
  function g_hs(model::sPCSAFTModel, V, T, z)
      η = @f(ζ,3)
      return (1-η/2)/(1-η)^3
  end
  
  function a_hs(model::sPCSAFTModel, V, T, z)
      η = @f(ζ,3)
      return (4η-3η^2)/(1-η)^2
  end
  ```

The rest works exactly as it would with the `PCSAFT` example.
