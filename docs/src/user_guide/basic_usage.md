# Clapeyron User Guide

Welcome to Clapeyron!

Once Clapeyron is installed, it can be loaded using:

```julia
using Clapeyron
```

We may create a model object by calling the constructor of the respective model. For example,

```julia
model1 = PCSAFT(["methanol", "ethanol"])
model2 = softSAFT(["ethane", "water"])
```

Clapeyron provides a database for a selection of parameters that are currently published in the literature. If you would like to use a custom database, you may link to it using the optional parameter `userlocations`.

```julia
model3 = PCSAFT(["my_comp1", "my_comp2"];
    userlocations = ["path/to/database_like.csv", "path/to/database_pair.csv"])
```

Clapeyron support three types of parameters: like, pair and assoc. More details can be found in the custom databases guide.

Once we have our model object, we will be able to call the respective thermodynamic methods to obtain the properties that we are looking for. For example, to find the isobaric heat capacity of a 0.5 mol methanol and 0.5 mol ethanol mixture using PC-SAFT at a pressure of 10 bar and a temperature of 300 K, we just call the `isobaric_heat_capacity(model, p, T, z)` function with the desired model and conditions as parameters.

```julia
Cp = isobaric_heat_capacity(model1, 10e5, 300, [0.5, 0.5])
```

The functions for the physical properties that we currently support are as follows:

```julia
volume(model, p, T, z)
sat_pure(model, T)
crit_pure(model)
enthalpy_vap(model, T)
pressure(model, V, T, z)
entropy(model, p, T, z)
chemical_potential(model, p, T, z)
internal_energy(model, p, T, z)
enthalpy(model, p, T, z)
Gibbs_free_energy(model, p, T, z)
Helmholtz_free_energy(model, p, T, z)
isochoric_heat_capacity(model, p, T, z)
isobaric_heat_capacity(model, p, T, z)
thermal_compressibility(model, p, T, z)
isentropic_compressibility(model, p, T, z)
speed_of_sound(model, p, T, z)
isobaric_expansitivity(model, p, T, z)
joule_thomson_coefficient(model, p, T, z)
compressibility_factor(model, p, T, z)
inversion_temperature(model, p, z)
```

`Clapeyron` also supports physical units through the use of `Unitful.jl`.

```julia
using Unitful
import Unitful: bar, °C, mol

Cp2 = isobaric_heat_capacity(model1, 5bar, 25°C, [0.5mol, 0.5mol])
```

Note that if you do not wish to import specific units, you may also just use a Unitful string, `pressure = 20u"psi"`.

We also support group-contribution models like SAFT-*ɣ* Mie. We have a database of species with the number of each group associated with it for easy lookup, but you may also use your own combinations. We use a tuple of the name of the molecule and an array of the group-multiplicity mappings.

```julia
model4 = SAFTgammaMie([
        "ethanol",
        ("nonadecanol", ["CH3"=>1, "CH2"=>18, "OH"=>1]),
        ("ibuprofen", ["CH3"=>3, "COOH"=>1, "aCCH"=>1, "aCCH2"=>1, "aCH"=>4])
        ]
    )
```

