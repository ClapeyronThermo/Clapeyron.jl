```@meta
CurrentModule = Clapeyron
```

Once Clapeyron is installed, it can be loaded using:

```julia
using Clapeyron
```

## Creating a thermodynamic model

We may create a model object by calling the constructor of the respective equation of state.
For example,

```julia
model1 = PCSAFT(["methanol"])
model2 = PR(["ethane", "water"])
model3 = GERG2008(["propane","pentane"])
```

## Group-contribution Models

We also support group-contribution models like SAFT‑$\gamma$ Mie.
We have a database of species with the number of each group associated with it for easy lookup, but you may also use your own combinations.
We use a tuple of the name of the molecule and an array of the group-multiplicity mappings.
For example

```julia
model4 = SAFTgammaMie([
        "ethanol",
        ("ibuprofen", ["CH3"=>3, "COOH"=>1, "aCCH"=>1, "aCCH2"=>1, "aCH"=>4])])
```

In some group-contribution approaches, one may need to specify some structural information (such as gc‑PC‑SAFT), such as the number of bonds between groups.
This can be done as follows:

```julia
model5 = gcPCSAFT([
        ("ethanol", ["CH3" => 1, "CH2OH" => 1], [("CH3", "CH2OH") => 1])
        ("octane", ["CH3" => 2, "CH2" => 6], [("CH3", "CH2") => 2, ("CH2", "CH2") => 5])])
```

## Available models

One can find out more about the information stored within these model objects in the API documentation.
In terms of equations of state available, we have the following default models:

**Cubics**:

- van der Waals ([`vdW`](@ref))
  - Clausius ([`Clausius`](@ref))
  - Berthelot ([`Berthelot`](@ref))
- Redlich–Kwong ([`RK`](@ref))
  - Soave–Redlich–Kwong ([`SRK`](@ref))
  - Predictive Soave–Redlich–Kwong ([`PSRK`](@ref))
  - Translated-and-Consistent Redlich–Kwong ([`tcRK`](@ref))
- Peng–Robinson ([`PR`](@ref))
  - Peng–Robinson (1978) ([`PR78`](@ref))
  - "Universal Mixing Rule" Peng–Robinson ([`UMRPR`](@ref))
  - Volume-Translated Peng–Robinson ([`VTPR`](@ref))
  - Translated-and-Consistent Peng–Robinson ([`tcPR`](@ref))
  - Consistent Peng–Robinson + Twu ([`cPR`](@ref))
  - Quantum Corrected Peng–Robinson ([`QCPR`](@ref))
  - Enhanced Predictive Peng–Robinson (1978) ([`EPPR78`](@ref))
- Patel–Teja ([`PatelTeja`](@ref))
  - Patel–Teja–Valderrama ([`PTV`](@ref))
- Kumar–Upadhyay ([`KU`](@ref))
- Redlich–Kwong–Peng–Robinson ([`RKPR`](@ref))

**SAFT**:

- SAFT ([`ogSAFT`](@ref))
- CK (Chen and Kreglewski) SAFT ([`CKSAFT`](@ref))
  - Simplified CK‑SAFT ([`sCKSAFT`](@ref))
- [`BACKSAFT`](@ref)
- Lennard-Jones SAFT ([`LJSAFT`](@ref))
- SAFT, Variable Range (VR), Square Well (SW) ([`SAFTVRSW`](@ref))
- Cubic plus Association ([`CPA`](@ref))
  - Simplified CPA ([`sCPA`](@ref))
- Soft SAFT, with Lennard-Jones function from Johnson et al. (1993) ([`softSAFT`](@ref))
- Soft SAFT, with Lennard-Jones function from Thol et al. (2016) ([`softSAFT2016`](@ref))
- Perturbed-Chain SAFT ([`PCSAFT`](@ref))
  - Simplified PC‑SAFT ([`sPCSAFT`](@ref))
  - PC‑SAFT with T-dependent kᵢⱼ and special correlation for water ([`pharmaPCSAFT`](@ref))
  - Heterogeneous GC‑PC‑SAFT ([`gcPCSAFT`](@ref))
  - PC‑SAFT with Gᴱ mixing rule ([`GEPCSAFT`](@ref))
- SAFT‑VR with Mie potential ([`SAFTVRMie`](@ref))
  - SAFT‑VR with quantum corrected Mie potential ([`SAFTVRQMie`](@ref))
- SAFT‑$\gamma$ Mie ([`SAFTgammaMie`](@ref))
  - Structural SAFT‑$\gamma$ Mie ([`structSAFTgammaMie`](@ref))

**Activity coefficient** (N.B. these models only provide VLE properties for mixtures):

- [`Wilson`](@ref)
- Non-random two-liquid ([`NRTL`](@ref))
  - NRTL, temperature-dependent interaction ([`aspenNRTL`](@ref))
- *Universal quasichemical Activity Coefficients* (UNIQUAC): ([`UNIQUAC`](@ref))
- *UNIQUAC Functional-group Activity Coefficients* (UNIFAC): ([`UNIFAC`](@ref))
  - UNIFAC‑FV ([`UNIFACFV`](@ref))
  - UNIFAC‑FV (polymer blends) ([`UNIFACFVPoly`](@ref))
- Conductor-like Screening Model Segment Activity Model (COSMO‑SAC)
  - COSMO‑SAC (2002 version) ([`COSMOSAC02`](@ref))
  - COSMO‑SAC (2010 version) ([`COSMOSAC10`](@ref))
  - COSMO‑SAC with dispersive interactions ([`COSMOSACdsp`](@ref))

**Empirical**:

- General MultiParameter Empiric model ([`MultiFluid`](@ref))
  - GERG‑2008 EoS for Natural Gas ([`GERG2008`](@ref))
  - EOS‑LNG for Liquefied Natural Gas ([`EOS_LNG`](@ref))
  - IAPWS‑95 Water reference ([`IAPWS95`](@ref))
  - Propane Reference ([`PropaneRef`](@ref))
  - Lennard-Jones Reference from Thol et al. (2016) ([`LJRef`](@ref))
  - Ammonia Reference (2023) ([`Ammonia2023`](@ref))
  - Multiparameter EoS + Activity ([`HelmAct`](@ref))

**Extended Corresponding States ([`ECS`](@ref))**:

- SPUNG ([`SPUNG`](@ref))

One can find out more about each of these equations of state within our background documentation.
Nevertheless, all of these equations are compatible with all methods available in our package.

There a few optional arguments available for these equations which will be explained below.
One of these is specifying the location of the parameter databases, the details of which can be found in our Custom databases documentation.

## Specifying an ideal term

Both SAFT and cubic-type equations of state rely upon an ideal model.
By default, Clapeyron uses what we refer to as the [`BasicIdeal`](@ref) model to account for the ideal contribution which does not require any parameters.
For properties which only have derivatives with respect to volume or composition (e.g. volume, isothermal compressibility, critical points, saturation points), or monoatomic species (e.g. noble gases), this is perfectly fine.
However, for any other properties or species, the results obtained will most likely be quite poor.
This is because this model does not account for the rotational and vibrational modes of the species.
To amend this, we provide three additional ideal models to be used instead:

- Monomer ideal correlation ([`MonomerIdeal`](@ref))
- Walker and Haslam's ideal correlation ([`WalkerIdeal`](@ref))
- Joback's ideal correlation ([`JobackIdeal`](@ref))
- Reid's polynomial correlation ([`ReidIdeal`](@ref))
- Aly-Lee's correlation ([`AlyLeeIdeal`](@ref))
- PPDS correlation ([`PPDSIdeal`](@ref))
- MultiParameter Empiric Ideal correlations ([`EmpiricIdeal`](@ref))

These can be specified for any of the SAFT or cubic-type equations of state using:

```julia
model5 = PCSAFT(["carbon dioxide"]; idealmodel = WalkerIdeal)
```

Everything else will work as normal.

## Specifying an alpha function

vdW, RK and PR cubic equations rely on an alpha function (SRK is technically just RK but with a different alpha function).
Whilst we use the defaults for both RK and PR, it is possible to toggle between them.
For example:

```julia
model6 = RK(["ethane","propane"];alpha=SoaveAlpha)
```

The above model would be equivalent to a model built by SRK directly.
We support the following alpha functions:

- [`RKAlpha`](@ref): This is the default alpha function for regular RK.
- [`SoaveAlpha`](@ref): This is the default alpha function for SRK.
- [`PRAlpha`](@ref): This is the default alpha function for regular PR.
- [`PR78Alpha`](@ref): This is the default alpha function for PR78.
- [`KUAlpha`](@ref): This is the default alpha function for KU
- [`RKPRAlpha`](@ref): This is the default alpha function for RKPR
- [`BMAlpha`](@ref): This is the modified alpha function proposed by Boston and Mathias designed to improve estimates above the critical point.
  This works for both PR and RK.
- [`TwuAlpha`](@ref): Proposed by Twu et al., this alpha function uses species-specific parameters rather than correlation and, thus, is slightly more accurate than regular alpha functions.
  It was intended to be used with PR and is used in VTPR, tcPR, cPR and tcRK.
- [`Twu88Alpha`](@ref): An earlier version of `TwuAlpha`, that uses 2 parameters instead of 3.
- [`MTAlpha`](@ref): Proposed by Magoulas and Tassios, this alpha function is essentially like the regular PR alpha function only to a higher order.
  It is used within UMRPR.
- [`Soave2019Alpha`](@ref): Updated soave correlations for `PR` and (s)`RK` models.

## Specifying a mixing rule

Only relevant to cubic equations of state and mixtures, we can alternate between different mixing rules in case these may result in better predictions.
We can toggle between these mixing rules:

```julia
model7 = RK(["ethane","propane"];mixing=KayRule)
```

We currently support:

- [`vdW1fRule`](@ref): The standard van der Waals one-fluid mixing rule which is the default in all cubics.
- [`KayRule`](@ref): Takes an approach closer to the mixing rules used in SAFT.
- [`HVRule`](@ref): The Huron–Vidal mixing rule with uses information from activity coefficient models to form the mixing rule.
  It is meant to be more accurate than regular mixing rules.
  As it requires an activity coefficient model, this must be specified:

  ```julia
  model7 = RK(["methanol","benzene"];mixing=HVRule,activity=Wilson)
  ```

- [`MHV1Rule`](@ref): The modified Huron–Vidal mixing rule proposed by Michelsen to first order.
  This has rather significant improvements over the regular mixing rule.
  Also needs an activity model to be specified.
- [`MHV2Rule`](@ref): The modified Huron–Vidal mixing rule proposed by Michelsen to second order.
  This is meant to be an improvement over the first order rule.
  Also needs an activity model to be specified.
- [`WSRule`](@ref): The Wong–Sandler mixing rule which also relies on an activity model.
  The equations are slightly more complicated but it is meant to be an improvement compared to `HVRule`.
  Also needs an activity model to be specified.
- [`modWSRule`](@ref): A modified Wong–Sandler mixing rule, that reduces to `vdW1fRule` when there is no nonideal mixtures.
- [`LCVMRule`](@ref): The Linear Combination of Vidal and Michelsen mixing rules is designed for asymmetric mixtures.
  Also needs an activity model to be specified.

If one goes looking within the source code, they will also find [`VTPRRule`](@ref), [`PSRKRule`](@ref),[`PPR78Rule`](@ref), [`QCPRRule`](@ref) and [`UMRRule`](@ref); these are only intended for use in their respective models and shouldn't be used otherwise.
However, it is still possible to toggle between them.

## Specifying a volume translation method

In order to improve the predictions of bulk properties in cubics, without affecting VLE properties, a volume translation method can be used which simply shifts the volume within the cubics by `c`.
The default for all cubics is `NoTranslation`, however, we can toggle between the methods:

```julia
model7 = RK(["ethane","propane"];translation=PenelouxTranslation)
```

We support the following methods:

- [`PenelouxTranslation`](@ref): Used in PSRK.
- [`RackettTranslation`](@ref): Used in VTPR.
- [`MTTranslation`](@ref): Used in UMRPR.
- [`ConstantTranslation`](@ref)

Note that not all these methods will be compatible with all species as they require the critical volume of the species.

## Using an Activity coefficient model

Activity coefficient models are primarily designed to obtain accurate estimate of mixture VLE properties *below* the critical point of all species.
Whilst not as flexible as other equations of state, they are computationally cheaper and, generally, more accurate.
The activity coefficients are obtained as only a function of temperature and composition ($\gamma (T,\mathbf{x})$), meaning we can simply use modified Raoult's law to obtain the bubble (and dew) point:

``y_ip= x_i\gamma_ip_{\mathrm{sat},i}``

The only problem here is that another model must provide the saturation pressure $p_{\mathrm{sat},i}$.
By default, this is chosen to be PR; however, one can toggle this setting as well:

```julia
model3 = UNIFAC(["methanol","benzene"];puremodel=PCSAFT)
```

Everything else will work as normal (so long as the species are also available within the specified pure model).

## Available properties

Once we have our model object, we will be able to call the respective thermodynamic methods to obtain the properties that we are looking for.
For example, to find the isobaric heat capacity of a 0.5 mol methanol and 0.5 mol ethanol mixture using PC‑SAFT at a pressure of 10 bar and a temperature of 300 K, we just call the `isobaric_heat_capacity(model, p, T, z)` function with the desired model and conditions as parameters.

```julia
Cp = isobaric_heat_capacity(model1, 10e5, 300, [0.5, 0.5])
```

The functions for the physical properties that we currently support are as follows:

- Bulk properties:

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
  alphaV = isobaric_expansivity(model, p, T, z)
  muJT = joule_thomson_coefficient(model, p, T, z)
  Z = compressibility_factor(model, p, T, z)
  gamma = activity_coefficients(model, p, T, z)
  ```

  All the above functions have two optional arguments (although, technically, z is an optional argument if you're only obtaining properties for a pure species):

  - `phase`: If you already know the phase of the species and want a (minor) speed-up, you can specify it.
    For example:

    ```julia
    V = volume(model, p, T, z; phase=:liquid)
    ```

    The default value is `:unknown` where it will find both the vapour and liquid roots first and determine which has the lowest Gibbs free energy.

  - `threaded`: This determines whether or not to run the vapour and liquid calculations in parallel or not and is only relevant for when the phases are unknown and non-cubic models.

    ```julia
    V = volume(model, p, T, z; threaded=false)
    ```

    The default value is `true`.
    This shouldn't change the results.

  Most of the above functions also accept the `vol0` optional keyword argument, which specifies an initial guess for the [Clapeyron.volume](@ref) solver.

  Note that all of the above functions can be broadcast i.e. if `T` is an array, instead of a for loop, we can simply:

  ```julia
  Cp = isobaric_heat_capacity.(model, p, T, z)
  ```

- Vapour–liquid, liquid–liquid and vapour–liquid–liquid equilibrium properties:

  - For pure species:

    ```julia
    (p_sat, V_l_sat, V_v_sat) = saturation_pressure(model, T)
    (T_sat, V_l_sat, V_v_sat) = saturation_temperature(model, p)
    H_vap = enthalpy_vap(model, T)
    ```

  - For mixtures:

    ```julia
    (p_bub, V_l_bub, V_v_bub, y) = bubble_pressure(model, T, x)
    (T_bub, V_l_bub, V_v_bub, y) = bubble_temperature(model, p, x)
    (p_dew, V_l_dew, V_v_dew, x) = dew_pressure(model, T, y)
    (T_dew, V_l_dew, V_v_dew, x) = dew_temperature(model, p, y)
    (p_LLE, V_l_LLE, V_ll_LLE, xx) = LLE_pressure(model, T, x)
    (T_LLE, V_l_LLE, V_ll_LLE, xx) = LLE_temperature(model, p, x)
    (p_az, V_l_az, V_v_az, x) = azeotrope_pressure(model, T)
    (T_az, V_l_az, V_v_az, x) = azeotrope_temperature(model, p)
    (p_VLLE,V_l_sat, V_ll_sat, V_v_sat, x, xx, y) = VLLE_pressure(model, T)
    (T_VLLE,V_l_sat, V_ll_sat, V_v_sat, x, xx, y) = VLLE_temperature(model, p)
    ```

  All the above arguments take in an optional argument for the initial guess:

  ```julia
  (p_sat, V_l_sat, V_v_sat) = saturation_pressure(model, T;v0=log10.([V_l0,V_v0]))
  ```

  Although our calculations tend to be quite robust, this argument is generally useful for when one wants to obtain smooth VLE envelopes quickly when making figures
  Here, you'd use a for loop where each iteration uses the previous' iteration value as an initial guess (except the first iteration).
  For example:

  ```julia
  (p_sat, V_l_sat, V_v_sat) = saturation_pressure(model, T[1])
  for i in 2:length(T)
    A = saturation_pressure(model,T[i];v0=log10.([V_l_sat[i-1],V_v_sat[i-1]]))
    append!(p_sat,A[1])
    append!(V_l_sat,A[2])
    append!(V_v_sat,A[3])
  end
  ```

- Critical properties:
  - For pure species:

    ```julia
    (T_c, p_c, V_c) = crit_pure(model)
    ```

  - For mixtures:

    ```julia
    (T_c, p_c, V_c) = crit_mix(model, z)
    (p_UCST, V_UCST, x_UCST) = UCST_mix(model, T)
    (T_UCEP, p_UCEP, V_l_UCEP, V_v_UCEP, x, y) = UCEP_mix(model)
    ```

  Like the above functions, for `crit_mix`, you can also specify initial guesses to produce smooth critical curves.
- Miscellaneous:

  ```julia
  T = inversion_temperature(model, p, z)
  B = second_virial_coefficient(model, T, z)
  ```

`Clapeyron` also supports physical units through the use of `Unitful.jl`.

```julia
using Unitful
import Unitful: bar, °C, mol, kg, l
model_unit = PCSAFT(["methanol","water"])
Cp2 = isobaric_heat_capacity(model_unit, 5bar, 25°C, [0.5kg, 0.5kg])  # isobaric heat capacity of 1 mol of mixture, at a pressure of 5 bar
Cp2 = isobaric_heat_capacity(model_unit, 1.0l/kg, 25°C, [0.4kg, 0.6kg])  # isobaric heat capacity of 1 kg of mixture, at a volume of 1 L/kg
```

Note that if you do not wish to import specific units, you may also just use a Unitful string, `pressure = 20u"psi"`.
This is only supported for bulk properties.
