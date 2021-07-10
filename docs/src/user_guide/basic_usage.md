# Clapeyron User Guide

Welcome to Clapeyron!

Once Clapeyron is installed, it can be loaded using:

```julia
using Clapeyron
```

We may create a model by calling the constructor of the respective model. For example,

```julia
model1 = PCSAFT(["methanol", "ethanol"])
model2 = SAFTgammaMie(["ethane", "water"])
```

Clapeyron provides a database for a selection of parameters that are currently published in the literature. If you would like to use a custom database, you may link to it using the optional parameter `userlocations`.

```julia
model = PCSAFT(["my_comp1", "my_comp2"];
    userlocations = ["path/to/database_like", "path/to/database_unlike"])
```



For example, to find the isobaric heat capacity at equal mole fractions at a pressure of 1 bar and a temperature of 300 K, we just write

```julia
Cp = get_isobaric_heat_capacity(model, create_z(model, [0.5, 0.5]), 10e5, 300)
```

The functions for the physical properties that we currently support are as follows:

```julia
get_volume(model, z, p, T)
get_sat_pure(model, T)
get_crit_pure(model)
get_enthalpy_vap(model, T)
get_pressure(model, z, v, T)
get_entropy(model, z, p, T)
get_chemical_potential(model, z, p, T)
get_internal_energy(model, z, p, T)
get_enthalpy(model, z, p, T)
get_Gibbs_free_energy(model, z, p, T)
get_Helmholtz_free_energy(model, z, p, T)
get_isochoric_heat_capacity(model, z, p, T)
get_isobaric_heat_capacity(model, z, p, T)
get_thermal_compressibility(model, z, p, T)
get_isentropic_compressibility(model, z, p, T)
get_speed_of_sound(model, z, p, T)
get_isobaric_expansitivity(model, z, p, T)
get_Joule_Thomson_coefficient(model, z, p, T)
```
