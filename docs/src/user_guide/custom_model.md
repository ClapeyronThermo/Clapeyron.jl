## Definitions

Within `Clapeyron`, it is possible to create your own model / equation of state whilst still using all of the property estimation tools we provide. In order to do this, you must create the model. As such, you need to define 3 things:

```julia
@newmodel ModelName ModelType ModelParam
```

1. Model name
2. Model type
3. Model parameters

We now give a brief overview of each of these aspects.

### Model name

This is simply the name of the equation of state which will be used to generate the model. It will be useful to use the same name for all other definitions within the model. Within `Clapeyron`, we've tried to keep these names as obvious as possible (_e.g._ `vdW`, `PR`, `PCSAFT`, `SAFTVRMie`, `GERG2008`, _etc._). For your own model, this can be whatever you want, but make sure to keep it consistent.

### Model type

Each model is assigned an abstract type (kind of like a label), the most general of which is `EoSModel`. From this parent, we branch into more-specific EoS types:

- `SAFTModel`: These are the models which have three parameters in common: segment size, $\sigma$, potential depth, $\epsilon$, and number of segments, $m$. All other SAFT-type models branch from this parent (`PCSAFTModel`, `SAFTVRMieModel`, `softSAFTModel`, _etc._). 
- `CubicModel`: These are the models whose parameters can be obtained from the critical temperature and pressure. With the exception of `CPAModel`, all cubics have a common structure where one can re-arrange the equation for the pressure as a third-order polynomial. As such, we define a subtype of `CubicModel`, `ABCubicModel` (_e.g._ `vdWModel`, `RKModel`, `SRKModel`, `PRModel`).
- `EmpiricHelmholtzModel`: These are the high-accuracy, multi-parameter models for specific species or systems (_e.g._ `GERG2008Model`, `IAPWS95Model`). There is no general structure to the models and they are treated as self-contained.
- `IdealModel`: Often overlooked, these models supplement the `SAFTModel` and `CubicModel` by providing the ideal contribution. Whilst the parameters and structure aren't usually the same between ideal models, this is unnecessary as the equation for the pressure is always $pV=Nk_\mathrm{B}T$ . 

A key aspect of making sure the methods available in `Clapeyron` will be compatible with your model is that the model must be a sub-type of either `SAFTModel` or `CubicModel`. This is because most of the methods (especially the initial guesses) are defined for a particular set of parameters or functional forms of the equations. If your model is not of either of these types, it will be necessary to define a few additional functions:

- `lb_volume(model::YourModel,T,z)`: This must output the smallest possible value for your particular model using the model parameters. In SAFT equations, this is equivalent to a packing fraction of one ($\eta=1$) whilst in cubics, it is equivalent to the `b` parameter.
- `T_scale(model::YourModel,z)`: This must output the temperature scaling for your model using the model parameters. In SAFT equations, this is usually the potential depth whilst in cubics it is the critical temperature.
- `p_scale(model::YourModel,z)`: This must output the pressure scaling for your model using the model parameters. In cubics, this is the critical pressure whilst in SAFT we use a more complicated definition using the segment sizes and potential depths.

The only exception to this are ideal models since these do not require numerical methods to be used and usually suppplement other models.

### Model parameters

The parameters for a particular system are all stored within a struct whose type is also linked to the model in the same way as above only we now use the `Param` suffix instead of `Model`. The hierarchy is also still the same. Each parameter within the struct is also assigned an abstract type. To explain this, we consider the generic parameter structs for SAFT-type and cubic-type models (for your own models, you can put whatever you need in these structs):

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

As we can see, we have a few new types to define for each of these parameters:

- `SingleParam`: These are the parameters associated with a pure species or indexed by a single index ($x_i$). For example, the molecular weight, number of segments, critical temperature and pressure.
- `PairParam`: These are the parameters associated with a pair of species or indexed by two indices ($x_{ij}$). If the two indices are the same ($i=j$), they refer to the like species, otherwise ($i\neq j$), they refer to the unlike species. For example, the cubic `a` and `b` parameters, the segment size (`sigma`) and potential depth (`epsilon`).
- `AssocParam`: There are the parameters associated with both a pair of species and association sites (see background documentation for what these are). They are indexed by four indices ($x_{ij,ab}$). For example, potential depth of the association interaction (`epsilon_assoc`) and the bonding volume (`bondvol`).

There are two additional parameter types needed for equations of state which are not stored within the EoSParam structs:

- `SiteParam`: Objects of this type will simply store the number of sites of a particular type on a species (_e.g._ how many H-sites on water?).
- `GroupParam`: Objects of this type are only relevant to the group contribution method of SAFT-$\gamma$ Mie. They will store the number of groups of a particular type within a species (_e.g._ propane has two CH$_3$ groups).

Unless you are doing something particularly new, you shouldn't need to interact with these functions.

## PC-SAFT Example

Once all the above has been defined, we are ready to build our own model. Let us imagine we are trying to implement `PCSAFT`. 

1. We first define the Model name, type and parameters:

```julia
# Defining the Model type
abstract type PCSAFTModel <: SAFTModel end

# Defining the Model params
struct PCSAFTParam <: EoSParam
    Mw::SingleParam{Float64}
    segment::SingleParam{Float64}
    sigma::PairParam{Float64}
    epsilon::PairParam{Float64}
    epsilon_assoc::AssocParam{Float64}
    bondvol::AssocParam{Float64}
end

# Creating the model
@newmodel PCSAFT PCSAFTModel PCSAFTParam

# Export the model
export PCSAFT
```

2. The next step is to write a function that builds our model:

```Julia
# All model functions should have the same header.
function PCSAFT(components::Array{String,1}; idealmodel=BasicIdeal, userlocations::Array{String,1}=String[], verbose=false)
  	# Obtain the parameters from the database. We can define some generic ones that are already available in Clapeyron. 
    params = Clapeyron.getparams(components, ["SAFT/PCSAFT","properties/molarmass.csv"]; userlocations=userlocations, verbose=verbose)
  
    # Obtain the various parameters from params.
    segment = params["m"]
    k = params["k"]
    Mw = params["Mw"]
    params["sigma"].values .*= 1E-10
  
    # In some cases, we may not have the unlike parameters and will need to use combining rules. You can also define your own combining rules for this.
    sigma = Clapeyron.sigma_LorentzBerthelot(params["sigma"])
    epsilon = Clapeyron.epsilon_LorentzBerthelot(params["epsilon"], k)
  
    epsilon_assoc = params["epsilon_assoc"]
    bondvol = params["bondvol"]
  
    # Build the sites object for associating species.
    sites = SiteParam(Dict("e" => params["n_e"], "H" => params["n_H"]))
  
    # Create the parameter structure.
    packagedparams = PCSAFTParam(Mw, segment, sigma, epsilon, epsilon_assoc, bondvol)
  
    # Although optional, it's generally good practise to cite your models!
    references = ["10.1021/ie0003887", "10.1021/ie010954d"]

    # Build the model.
    model = PCSAFT(packagedparams, sites, idealmodel; references=references, verbose=verbose)
  
    # Output.
    return model
end
```

3. Define all of the model equations. All functions will need to have the `Clapeyron.` suffix. We recommend that you use UNICODE characters wherever possible!

   ```julia
   function Clapeyron.a_res(model::PCSAFTModel, V, T, z)
       return @f(a_hc) + @f(a_disp) + @f(a_assoc)
   end
   
   function Clapeyron.a_hc(model::PCSAFTModel, V, T, z)
       x = z/∑(z)
       m = model.params.segment.values
       m̄ = ∑(x .* m)
       return m̄*@f(a_hs) - ∑(x[i]*(m[i]-1)*log(@f(g_hs,i,i)) for i ∈ @comps)
   end
   
   function Clapeyron.d(model::PCSAFTModel, V, T, z, i)
       ϵii = model.params.epsilon.diagvalues[i]
       σii = model.params.sigma.diagvalues[i]
       return σii * (1 - 0.12exp(-3ϵii/T))
   end
   
   function Clapeyron.ζ(model::PCSAFTModel, V, T, z, n)
       ∑z = ∑(z)
       x = z * (one(∑z)/∑z)
       m = model.params.segment.values
       res = N_A*∑z*π/6/V * ∑((x[i]*m[i]*@f(d,i)^n for i ∈ @comps))
   end
   
   function Clapeyron.g_hs(model::PCSAFTModel, V, T, z, i, j)
       di = @f(d,i)
       dj = @f(d,j)
       ζ2 = @f(ζ,2)
       ζ3 = @f(ζ,3)
       return 1/(1-ζ3) + di*dj/(di+dj)*3ζ2/(1-ζ3)^2 + (di*dj/(di+dj))^2*2ζ2^2/(1-ζ3)^3
   end
   
   function Clapeyron.a_hs(model::PCSAFTModel, V, T, z)
       ζ0 = @f(ζ,0)
       ζ1 = @f(ζ,1)
       ζ2 = @f(ζ,2)
       ζ3 = @f(ζ,3)
       return 1/ζ0 * (3ζ1*ζ2/(1-ζ3) + ζ2^3/(ζ3*(1-ζ3)^2) + (ζ2^3/ζ3^2-ζ0)*log(1-ζ3))
   end
   
   # INSERT REST OF CODE
   ```

   You will notice in the above the macro `@f`. This macro will evaluate functions which have the inputs `fun(model, V, T, z)`; any additional inputs can also be added. 

4. With all the above defined in a single script, we can save the file as `PCSAFT.jl` and then include it in our jupyter notebooks (for example) and use the model with all of our existing method:

   ```julia
   include("PCSAFT.jl")
   
   model = PCSAFT(["carbon dioxide"])
   
   p = 20e6
   T = range(290,460,length=200)
   
   Cp = isobaric_heat_capacity.(model, p, T)
   
   (T_c, p_c, V_c) = crit_pure(model)
   
   T_sat = range(220,T_c,length=200)
   
   (p_sat, V_l_sat, V_v_sat) = sat_pure(model,T_sat)
   ```

   

