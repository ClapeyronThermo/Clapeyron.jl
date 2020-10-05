# OpenSAFT User Guide

Welcome to OpenSAFT!

Once OpenSAFT is installed, load it using

```julia
using OpenSAFT
```

We may create a model using the ```system(components, method; kwargs)``` function.

```julia
model = system(["methanol", "ethanol"], "PCSAFT")
```

OpenSAFT provides a database for a lot of paramaters that are currently publised in the literature. If you would like to use a custom database, copy the relevant files (in this case, PCSAFT) from the ```database``` directory in the package root directory to use as a template, and place it wherever you want, then create the model using

```julia
model = system(["methanol", "ethanol"], "PCSAFT";
    customdatabase_like = "{PATH}",
    customdatabase_unlike = "{PATH}",
    customdatabase_assoc = "{PATH}")
```

You may now find physical properties for this system!

Note that the composition has to be defined using the NamedArray type. We currently have a function called ```create_z(model, composition)``` to do this for you, but this will be converted into a macro later.

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
