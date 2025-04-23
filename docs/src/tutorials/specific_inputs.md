# Specific User Inputs

As seen before, the default inputs to compute any bulk property using any `EoSModel` are pressure, temperature, and moles. However, there is also an option to use alternative inputs. The specific combinantions are as follows:

1. Pressure and enthalpy `ph_flash`
2. Pressure and entropy `ps_flash`
3. Vapour fraction and pressure `qp_flash`
4. Vapour fraction and temperature `qt_flash`
5. Temperature and entropy `ts_flash`
6. Volume and temperature `vt_flash`

The following example demonstrates the use of `ph_flash`, but the same procedure applies to all of the user-specified input combinations above:

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


