```@meta
CurrentModule = Clapeyron
```

## Index

```@index
Pages = ["activity.md"]
```

# Activity Models

There are two alternatives on the definition of an activity model:
- Defining an excess Gibbs energy function
- Defining an activity coefficient function

those two can be converted between one form to another via:

``\gamma_{i} = \frac{\partial{G^E}}{\partial{n_i}}``

``\ln{\gamma_{i}} = \frac{1}{RT}\frac{\partial{G^E}}{\partial{n_i}}``

When defining one form, the other is derived automatically.

Those functions can also be derived from any arbitrary equation of state:

``\frac{\partial{G^E}}{\partial{n_i}}= \mu_i - \mu^{0}_i``

Where ``\mu_i`` and ``\mu^{0}_i`` are the mixture and pure chemical potentials of component ``i``. in this case, those potentials are dependent of the pressure. whereas activity models are usually only temperature and composition dependent.

```@docs
Clapeyron.Wilson
Clapeyron.NRTL
Clapeyron.aspenNRTL
Clapeyron.UNIQUAC
Clapeyron.ogUNIFAC
Clapeyron.UNIFAC
Clapeyron.PSRKUNIFAC
Clapeyron.VTPRUNIFAC
Clapeyron.UNIFACFV
Clapeyron.UNIFACFVPoly
Clapeyron.tcPRWilsonRes
Clapeyron.COSMOSAC02
Clapeyron.COSMOSAC10
Clapeyron.COSMOSACdsp
```
