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

We also support group-contribution models like SAFT-*ɣ* Mie. We have a database of species with the number of each group associated with it for easy lookup, but you may also use your own combinations. We use a tuple of the name of the molecule and an array of the group-multiplicity mappings.

```julia
model3 = SAFTgammaMie([
        "ethanol",
        ("nonadecanol", ["CH3"=>1, "CH2"=>18, "OH"=>1]),
        ("ibuprofen", ["CH3"=>3, "COOH"=>1, "aCCH"=>1, "aCCH2"=>1, "aCH"=>4])
        ]
    )
```

### Available equations of state

Cubic-type:

- van der Waals (`vdW`)
- Redlich-Kwong (`RK`)
- Soave-Redlich-Kwong (`SRK`)
- Peng-Robinson (`PR`)
- Cubic-plus-association (`CPA`)

SAFT-type:

- PC-SAFT (`PCSAFT`)
  - sPC-SAFT (`sPCSAFT`)
- SAFT-*ɣ* Mie (`SAFTgammaMie`)
- SAFT-VR Mie (`SAFTVRMie`)
  - SAFT-VRQ Mie (`SAFTVRQMie`)
- soft-SAFT (`softSAFT`)
- Original SAFT (`ogSAFT`)
- CK-SAFT (`CKSAFT`)
  - sCK-SAFT (`sCKSAFT`)
- LJ-SAFT (`LJSAFT`)
- BACK-SAFT (`BACKSAFT`)
- SAFT-VR SW (`SAFTVRSW`)

Multi-parameter equations:

- IAPWS-95 (`IAWPS95`)
- GERG-2008 (`GERG2008`)
- Propane Reference (`PropaneRef`)

We also support the SPUNG method (`SPUNG`)

### Specifying an ideal term

By default, Clapeyron uses what we refer to as the `BasicIdeal` model to account for the ideal contribution. For properties which only have derivatives with respect to volume or composition (_e.g._ pressure, isothermal compressibility, critical points, saturation points), or monoatomic species (_e.g._ noble gases), this is perfectly fine. However, for any other properties or species, the results obtained will most likely be quite poor. This is because this model does not account for the rotational and vibrational modes of the species. To amend this, we provide two additional ideal models to be used instead (more soon to come):

- Walker and Haslam's ideal correlation (`WalkerIdeal`)
- Reid's polynomial correlation (`ReidIdeal`)

These can be specified for any of the cubic or SAFT-type equations of state using:

```julia
model4 = PCSAFT(["carbon dioxide"]; idealmodel = WalkerIdeal)
```

Everything else will work as normal.

### Available properties

Once we have our model object, we will be able to call the respective thermodynamic methods to obtain the properties that we are looking for. For example, to find the isobaric heat capacity of a 0.5 mol methanol and 0.5 mol ethanol mixture using PC-SAFT at a pressure of 10 bar and a temperature of 300 K, we just call the `isobaric_heat_capacity(model, p, T, z)` function with the desired model and conditions as parameters.

```julia
Cp = isobaric_heat_capacity(model1, 10e5, 300, [0.5, 0.5])
```

The functions for the physical properties that we currently support are as follows:

- Bulk properties (must specify _p, T, z_):

  ```julia
  V = volume(model, p, T, z)
  p = pressure(model, V, T, z)
  S = entropy(model, p, T, z)
  mu = chemical_potential(model, p, T, z)
  U = internal_energy(model, p, T, z)
  H = enthalpy(model, p, T, z)
  G = Gibbs_free_energy(model, p, T, z)
  A = Helmholtz_free_energy(model, p, T, z)
  Cv = isochoric_heat_capacity(model, p, T, z)
  Cp = isobaric_heat_capacity(model, p, T, z)
  betaT = thermal_compressibility(model, p, T, z)
  betaS = isentropic_compressibility(model, p, T, z)
  u = speed_of_sound(model, p, T, z)
  alphaV = isobaric_expansitivity(model, p, T, z)
  muJT = joule_thomson_coefficient(model, p, T, z)
  Z = compressibility_factor(model, p, T, z)
  ```

  Note that all of the above functions can be broadcast _i.e._ if `T` is an array, instead of a for loop, we can simply:

  ```julia
  Cp = isobaric_heat_capacity.(model, p, T, z)
  ```

- Vapour-liquid equilibrium properties (must specify _T, z_):

  ```julia
  (p_sat, V_l_sat, V_v_sat) = sat_pure(model, T)
  H_vap = enthalpy_vap(model, T)
  ```

  Mixture VLE properties are currently in development

- Critical properties:

  ```julia
  (T_c, p_c, V_c) = crit_pure(model)
  ```

- Miscellaneous:

  ```julia
  T = inversion_temperature(model, p, z)
  B = second_virial_coefficient(model, T, z)
  ```

We note that whilst the bulk properties are all compatible with mixtures, the VLE and critical properties methods are only compatible with pure systems. The extension to mixtures is currently in development.

`Clapeyron` also supports physical units through the use of `Unitful.jl`.

```julia
using Unitful
import Unitful: bar, °C, mol

Cp2 = isobaric_heat_capacity(model1, 5bar, 25°C, [0.5mol, 0.5mol])
```

Note that if you do not wish to import specific units, you may also just use a Unitful string, `pressure = 20u"psi"`.

## Customisation

Although we provide many models, methods and parameters, Clapeyron also allows for easy customisation in all three of these aspects. To find out more how to customise your models, please read the relevant sections in the documentation.
