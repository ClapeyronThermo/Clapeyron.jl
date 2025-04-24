# Other Flash Computations

As seen before, the default inputs to compute any bulk property using any `EoSModel` are pressure, temperature, and moles. However, there is also an option to use alternative inputs. The specific combinantions are as follows:

1. Pressure and enthalpy `ph_flash`
2. Pressure and entropy `ps_flash`
3. Vapour fraction and pressure `qp_flash`
4. Vapour fraction and temperature `qt_flash`
5. Temperature and entropy `ts_flash`
6. Volume and temperature `vt_flash`


## Using P-H flash
The following example demonstrates the use of `ph_flash`, but the same procedure applies to all flash functions:
```julia
julia> model = cPR(["ethane","methane"],idealmodel = ReidIdeal);
julia> z = [1.0,1.0]; p = 101325; h = 100;
julia> flash_result = ph_flash(model,p,h,z)
Flash result at T = 299.938, p = 101325.0 with 1 phase:
 (x = [0.5, 0.5], β = 2.0, v = 0.0244958)
```
Once the flash_result is computed, other bulk properties can be determined as follows:

```julia
julia> s = entropy(model,flash_result)
-66.39869200962218

julia> mass_density(model,flash_result)
0.941244521997331
```
Additionally there are convinent modules which can be used to by-pass the manual computation of the flash result. Example: 

```julia
julia> using Clapeyron

julia> import Clapeyron: PH

julia> model = cPR(["ethane","methane"],idealmodel = ReidIdeal);

julia> z = [1.0,1.0]; p = 101325; h = 100;

julia> PH.entropy(model,p,h,z)
-66.39869200962218
```

Currently, convenience modules exist for PH, PS and VT flashes.

# Examples of using other flash computations

Below is an example of each flash computation for a 1:1 molar mixture of isopentane and isobutane:

```julia
julia> model = cPR(["isopentane","isobutane"],idealmodel = ReidIdeal)                                                    
PR{ReidIdeal, TwuAlpha, NoTranslation, vdW1fRule} with 2 components:
 "isopentane"
 "isobutane"
Contains parameters: a, b, Tc, Pc, Mw
```

## P-S flash
Using the `PS` module:
```julia
julia> z = [1.0,1.0]; p = 101325; s = 100;

julia> import Clapeyron: PS

julia> PS.temperature(model,p,s,z)
541.1993196556604

julia> PS.enthalpy(model,p,s,z)
69569.1104222583
```

## Q-P flash
Compute the entropy at vapour fraction 0.5 and pressure 101 325 Pa:
```julia
julia> z = [1.0,1.0]; p = 101325; q = 0.5;

julia> flash_result = qp_flash(model,q,p,z)
Flash result at T = 280.803, p = 101325.0 with 2 phases:
 (x = [0.667227, 0.332773], β = 1.0, v = 0.000105263)
 (x = [0.332773, 0.667227], β = 1.0, v = 0.0222189)

julia> entropy(model,flash_result)
-164.74025465755165
```

## Q-T flash
Entropy at vapour fraction 0.5 and temperature 300 K:

```julia
julia> z = [1.0,1.0]; T = 300; q = 0.5;

julia> flash_result = qt_flash(model,q,T,z)
Flash result at T = 300.0, p = 1.94999e5 with 2 phases:
 (x = [0.649352, 0.350648], β = 1.0, v = 0.000108864)
 (x = [0.350648, 0.649352], β = 1.0, v = 0.0120351)

julia> entropy(model,flash_result)
-153.12015827330828

julia> pressure(model,flash_result)
194998.54983747654
```

## T-S flash
Enthalpy and pressure at entropy –215 J/K and temperature 310 K: 

```julia
julia> z = [1.0,1.0]; T = 310; s = -215;

julia> flash_result = ts_flash(model,T,s,z)
Flash result at T = 310.0, p = 1.87265e6 with 1 phase:
 (x = [0.5, 0.5], β = 2.0, v = 0.000108424)

 julia> enthalpy(model,flash_result)
-41092.06962844136

julia> pressure(model,flash_result)
1.8726539417569228e6
```
## V-T flash
Enthalpy and pressure at volume 0.04 m³ and temperature 300 K:

```julia
julia> z = [1.0,1.0]; T = 300; v = 0.04;

julia> flash_result = vt_flash(model,v,T,z)
Flash result at T = 300.0, p = 1.19954e5 with 1 phase:
 (x = [0.5, 0.5], β = 2.0, v = 0.02)

julia> pressure(model,flash_result)
119953.80563632645

julia> enthalpy(model,flash_result)
-78.48634320658675
```
You can also use the VT module directly:
 ```julia
julia> import Clapeyron: VT

julia> VT.enthalpy(model,v,T,z)
-78.48634320658675
 ```



