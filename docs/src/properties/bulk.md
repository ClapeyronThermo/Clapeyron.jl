```@meta
CurrentModule = Clapeyron
```

## Contents

```@contents
Pages = ["bulk.md"]
```

## Index

```@index
Pages = ["bulk.md"]
```

## Volume–Temperature Based Properties

```@docs
Clapeyron.pressure
Clapeyron.second_virial_coefficient
Clapeyron.pip
```

## Pressure–Temperature Based Bulk Properties

In general almost all bulk properties follow the pattern:

```julia
function property(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)
    V = volume(model, p, T, z; phase, threaded, vol0)
    return VT_property(model,V,T,z)
end
```

So, you can calculate the property with Volume–Temperature variables by calling `VT_property(model,V,T,z)`.
Another way to do this is by using units, provided by `Unitful.jl`:

```julia
using Unitful
r = 18u"kg/m^3"
T = 373.15"K"
prop = helmholtz_free_energy(model,r,T,z,output = u"kJ")
```

Where `r` could be any molar or mass density, molar or mass volume, total volume or pressure.
It also supports mass and mol amounts defined as units for the composition (`z`).
If no units are provided for the composition, they will be considered moles.

### Methods that require first order VT derivatives

```@docs
Clapeyron.volume
Clapeyron.helmholtz_free_energy
Clapeyron.helmholtz_free_energy_res
Clapeyron.molar_density
Clapeyron.mass_density
Clapeyron.compressibility_factor
Clapeyron.gibbs_free_energy
Clapeyron.gibbs_free_energy_res
Clapeyron.entropy
Clapeyron.entropy_res
Clapeyron.enthalpy
Clapeyron.internal_energy
Clapeyron.internal_energy_res
```

### Methods that require second order VT derivatives

```@docs
Clapeyron.isochoric_heat_capacity
Clapeyron.isobaric_heat_capacity
Clapeyron.adiabatic_index
Clapeyron.isothermal_compressibility
Clapeyron.isentropic_compressibility
Clapeyron.speed_of_sound
Clapeyron.isobaric_expansivity
Clapeyron.joule_thomson_coefficient
```

### Methods that require first order composition derivatives

```@docs
Clapeyron.chemical_potential
Clapeyron.chemical_potential_res
Clapeyron.fugacity_coefficient
```

### Activity Coefficient

```@docs
Clapeyron.reference_chemical_potential
Clapeyron.reference_chemical_potential_type
Clapeyron.activity_coefficient
Clapeyron.activity
Clapeyron.aqueous_activity
```

### Mixing

```@docs
Clapeyron.mixing
```

## Initial guess functions

These methods are considered internal, they don't support `Symbolics.jl` or `Unitful.jl` overloads.

```@docs
Clapeyron.lb_volume
Clapeyron.T_scale
Clapeyron.p_scale
Clapeyron.x0_volume
Clapeyron.x0_volume_solid
Clapeyron.x0_volume_liquid
Clapeyron.x0_volume_gas
Clapeyron.volume_virial
Clapeyron.x0_sat_pure
Clapeyron.x0_psat
Clapeyron.x0_saturation_temperature
Clapeyron.antoine_coef
Clapeyron.x0_crit_pure
```
