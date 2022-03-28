```@meta
CurrentModule = Clapeyron
```

## Contents

```@contents
Pages = ["properties.md"]
```

## Index

```@index
Pages = ["properties.md"]
```

## Volume-Temperature Based Properties

```@docs
Clapeyron.pressure
Clapeyron.second_virial_coefficient
Clapeyron.pip
```

## Pressure-Temperature Based Bulk Properties

In general almost all bulk properties the pattern:
```julia
function property(model::EoSModel, p, T, z=SA[1.]; phase = :unknown,threaded=true)
    V = volume(model, p, T, z; phase=phase, threaded=threaded)
    return VT_property(model,V,T,z)
end
```
So, you can call the Volume-Temperature based version just by calling `VT_property(model,V,T,z).`
Another way to do this is by using units by `Unitful.jl`:
```julia
using Unitful
r = 18u"kg/m^3"
T = 373.15"K"
prop = helholtz_free_energy(model,r,T,z,output = u"kJ")
```
Where `r` could be any molar or mass density, molar or mass volume, total volume or pressure. it also supports mass and mol amounts defined as units in the `z` variable. if no units are provided by `z`, it will be assumed that they are moles.

### Methods that require first order VT derivatives
```@docs
Clapeyron.volume
Clapeyron.helmholtz_free_energy
Clapeyron.molar_density
Clapeyron.mass_density
Clapeyron.compressibility_factor
Clapeyron.gibbs_free_energy
Clapeyron.entropy
Clapeyron.entropy_res
Clapeyron.enthalpy
Clapeyron.internal_energy
```

### Methods that require second order VT derivatives
```@docs
Clapeyron.isochoric_heat_capacity
Clapeyron.isobaric_heat_capacity
Clapeyron.isothermal_compressibility
Clapeyron.isentropic_compressibility
Clapeyron.speed_of_sound
Clapeyron.isobaric_expansivity
Clapeyron.joule_thomson_coefficient
Clapeyron.compressibility_factor
```

### Methods that first order composition derivatives
```@docs
Clapeyron.chemical_potential
Clapeyron.chemical_potential_res
Clapeyron.fugacity_coefficient
```

### Misc
```
Clapeyron.mixing
```

## Single component properties

```@docs
Clapeyron.saturation_pressure
Clapeyron.enthalpy_vap
Clapeyron.crit_pure
Clapeyron.acentric_factor
```

## Multi component properties

```@docs
Clapeyron.bubble_pressure
Clapeyron.bubble_temperature
Clapeyron.dew_pressure
Clapeyron.dew_temperature
Clapeyron.azeotrope_pressure
Clapeyron.azeotrope_temperature
Clapeyron.LLE_pressure
Clapeyron.LLE_temperature
Clapeyron.VLLE_pressure
Clapeyron.VLLE_temperature
Clapeyron.crit_mix
Clapeyron.UCEP_mix
Clapeyron.UCST_mix
Clapeyron.gibbs_solvation
```

## Consistency and Stability

```@docs
Clapeyron.gibbs_duhem
Clapeyron.isstable
Clapeyron.mechanical_stability
Clapeyron.diffusive_stability
```

## TP Flash

```@docs
Clapeyron.tp_flash
Clapeyron.TPFlashMethod
Clapeyron.DETPFlash
Clapeyron.RRTPFlash
Clapeyron.numphases
```

## Initial guess functions

```@docs
Clapeyron.lb_volume
Clapeyron.T_scale
Clapeyron.p_scale
Clapeyron.x0_volume
Clapeyron.x0_volume_liquid
Clapeyron.x0_volume_gas
Clapeyron.volume_virial
Clapeyron.x0_sat_pure
Clapeyron.x0_crit_pure
```



