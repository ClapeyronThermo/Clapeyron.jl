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

## Bulk properties

By default, `Clapeyron.jl` exports functions that can be evaluated in a direct pressure-temperature (P-T-0) basis (also available in the `PT0` sub module.)

The basis is direct in the sense that only the most stable volume phase is calculated at the pressure-temperature-composition state, but no attempt is done to check if one or more phases can be formed at the input conditions.

A summary table can be found below:

| Property name                          | Main function                                   | Variants (full function name with descriptor)                                                                 |
|---------------------------------------|------------------------------------------------|----------------------------------------------------------------------------------------------------------------|
| Volume (total)                        | [`volume`](@ref)                               | –                                                                                                              |
| Temperature                            | [`temperature`](@ref)                             | –                                                                                                              |
| Pressure                               | [`pressure`](@ref)                             | –                                                                                                              |
| Helmholtz energy (total)               | [`helmholtz_energy`](@ref)                     | [`helmholtz_energy_res`](@ref) (residual), [`mass_helmholtz_energy`](@ref) (mass-based), [`helmholtz_free_energy`](@ref) (deprecated), [`helmholtz_free_energy_res`](@ref) (deprecated), [`mass_helmholtz_free_energy`](@ref) (deprecated) |
| Molar density                          | [`molar_density`](@ref)                        | –                                                                                                              |
| Mass density                           | [`mass_density`](@ref)                         | –                                                                                                              |
| Compressibility factor                 | [`compressibility_factor`](@ref)               | –                                                                                                              |
| Gibbs energy (total)                   | [`gibbs_energy`](@ref)                         | [`gibbs_energy_res`](@ref) (residual), [`mass_gibbs_energy`](@ref) (mass-based), [`gibbs_free_energy`](@ref) (deprecated), [`gibbs_free_energy_res`](@ref) (deprecated), [`mass_gibbs_free_energy`](@ref) (deprecated) |
| Entropy (total)                        | [`entropy`](@ref)                              | [`entropy_res`](@ref) (residual), [`mass_entropy`](@ref) (mass-based)                                         |
| Enthalpy (total)                       | [`enthalpy`](@ref)                             | [`enthalpy_res`](@ref) (residual), [`mass_enthalpy`](@ref) (mass-based)                                       |
| Internal energy (total)                | [`internal_energy`](@ref)                      | [`internal_energy_res`](@ref) (residual), [`mass_internal_energy`](@ref) (mass-based)                         |
| Isochoric heat capacity (molar)        | [`isochoric_heat_capacity`](@ref)              | [`mass_isochoric_heat_capacity`](@ref) (mass-based)                                                           |
| Isobaric heat capacity (molar)         | [`isobaric_heat_capacity`](@ref)               | [`mass_isobaric_heat_capacity`](@ref) (mass-based)                                                            |
| Adiabatic index                        | [`adiabatic_index`](@ref)                      | –                                                                                                              |
| Isothermal compressibility             | [`isothermal_compressibility`](@ref)           | –                                                                                                              |
| Isentropic compressibility             | [`isentropic_compressibility`](@ref)           | –                                                                                                              |
| Speed of sound                         | [`speed_of_sound`](@ref)                       | –                                                                                                              |
| Isobaric expansivity                   | [`isobaric_expansivity`](@ref)                 | –                                                                                                              |
| Joule‑Thomson coefficient              | [`joule_thomson_coefficient`](@ref)            | –                                                                                                              |
| Inversion temperature                  | [`inversion_temperature`](@ref)                | –                                                                                                              |
| Fundamental derivative of gas dynamics | [`fundamental_derivative_of_gas_dynamics`](@ref)| –                                                                                                              |
| Second virial coefficient              | [`second_virial_coefficient`](@ref)            | –                                                                                                              |
| Cross second virial coefficient        | [`cross_second_virial`](@ref)                  | –                                                                                                              |
| Equivolumetric cross second virial     | [`equivol_cross_second_virial`](@ref)          | –                                                                                                              |
| Chemical potential                     | [`chemical_potential`](@ref)                   | [`chemical_potential_res`](@ref) (residual)                                                                   |
| Fugacity coefficient                   | [`fugacity_coefficient`](@ref)                 | –                                                                                                              |
| Activity coefficient                   | [`activity_coefficient`](@ref)                 | –                                                                                                              |
| Activity                               | [`activity`](@ref)                             | –                                                                                                              |
| Aqueous activity                       | [`aqueous_activity`](@ref)                     | –                                                                                                              |
| Reference chemical potential           | [`reference_chemical_potential`](@ref)         | –                                                                                                              |
| Reference chemical potential type      | [`reference_chemical_potential_type`](@ref)    | –                                                                                                              |
| Mixing property                        | [`mixing`](@ref)                               | –                                                                                                              |
| Excess property                        | [`excess`](@ref)                               | –                                                                                                              |
| Partial property                       | [`partial_property`](@ref)                     | –                                                                                                              |
| Thermodynamic factor                   | [`thermodynamic_factor`](@ref)                 | –                                                                                                              |
| Phase identification parameter         | [`pip`](@ref)                                  | –                                                                                                              |
| Phase identification                   | [`identify_phase`](@ref)                       | –                                                                                                              |

### First order bulk methods

The following methods are calculated as combinations of first order derivatives of the helmholtz energy with respect to volume and temperature, or first order derivatives of the gibbs energy with respect to pressure and temperature.

```@docs
Clapeyron.volume
Clapeyron.pressure
Clapeyron.helmholtz_energy
Clapeyron.helmholtz_energy_res
Clapeyron.molar_density
Clapeyron.mass_density
Clapeyron.compressibility_factor
Clapeyron.gibbs_energy
Clapeyron.gibbs_energy_res
Clapeyron.entropy
Clapeyron.entropy_res
Clapeyron.enthalpy
Clapeyron.enthalpy_res
Clapeyron.internal_energy
Clapeyron.internal_energy_res
```

## Mass-based first order bulk properties

```@docs
Clapeyron.mass_enthalpy
Clapeyron.mass_entropy
Clapeyron.mass_internal_energy
Clapeyron.mass_gibbs_energy
Clapeyron.mass_helmholtz_energy
```

### Second order bulk methods

The following methods are calculated as combinations all derivatives up to the second order of the helmholtz energy with respect to volume and temperature, or all derivatives up to the second order of the gibbs energy with respect to pressure and temperature.

```@docs
Clapeyron.isochoric_heat_capacity
Clapeyron.isobaric_heat_capacity
Clapeyron.adiabatic_index
Clapeyron.isothermal_compressibility
Clapeyron.isentropic_compressibility
Clapeyron.speed_of_sound
Clapeyron.isobaric_expansivity
Clapeyron.joule_thomson_coefficient
Clapeyron.inversion_temperature
Clapeyron.fundamental_derivative_of_gas_dynamics
```

## Mass-based second order bulk properties

```@docs
Clapeyron.mass_isochoric_heat_capacity
Clapeyron.mass_isobaric_heat_capacity
```

## Temperature-only based properties

The following methods depend only on temperature and/or composition.

```@docs
Clapeyron.second_virial_coefficient
Clapeyron.cross_second_virial
Clapeyron.equivol_cross_second_virial
```

### Chemical potential functions

```@docs
Clapeyron.chemical_potential
Clapeyron.chemical_potential_res
Clapeyron.fugacity_coefficient
Clapeyron.activity_coefficient
Clapeyron.activity
Clapeyron.aqueous_activity
Clapeyron.reference_chemical_potential
Clapeyron.reference_chemical_potential_type
```

### Mixing and partial properties

```@docs
Clapeyron.mixing
Clapeyron.excess
Clapeyron.partial_property
Clapeyron.shape_factors
Clapeyron.thermodynamic_factor
```

### Phase identification

```@docs
Clapeyron.pip
Clapeyron.identify_phase
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
Clapeyron.x0_crit_pure
```

## Inverse property solvers

```@docs
Clapeyron.Tproperty
Clapeyron.Pproperty
```

### Bulk properties in other basis

```@docs
Clapeyron.VT0
Clapeyron.PT0
Clapeyron.PT
Clapeyron.VT
Clapeyron.PS
Clapeyron.PH
Clapeyron.QT
Clapeyron.QP
```


## Old aliases

The following functions were renamed to get rid of the `_free` naming. They still can be used, but they are considered deprecated.

```@docs
Clapeyron.gibbs_free_energy
Clapeyron.helmholtz_free_energy
Clapeyron.gibbs_free_energy_res
Clapeyron.helmholtz_free_energy_res
Clapeyron.mass_gibbs_free_energy
Clapeyron.mass_helmholtz_free_energy
```

