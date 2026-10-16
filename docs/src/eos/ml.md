# ML-Based Models

Most of the following models (with the exception of `UNIFAC2` models) are provided by [`MLThermoProperties.jl`](https://github.com/se-schmitt/MLThermoProperties.jl) and require the package to be installed.

## mod. UNIFAC 2.0 and UNIFAC 2.0

UNIFAC 2.0 and mod. UNIFAC 2.0 are enhanced versions of the classical group-contribution methods UNIFAC and mod. UNIFAC (Dortmund), respectively.
Missing interaction parameters are predicted using matrix completion, which significantly extends the applicability of the methods and leads to a higher prediction accuracy compared to the original versions.

```@docs ;canonical = false
Clapeyron.UNIFAC2
Clapeyron.ogUNIFAC2
```

## HANNA activity models

`HANNA` is a hard-constraint neural network model for the excess Gibbs energy `g^E` that predicts activity coefficients in a strictly thermodynamically consistent manner.
It only requires the SMILES of the components and the temperature as input.
The model satisfies thermodynamic boundary conditions by construction, ensuring consistency of the predicted activity coefficients.

Two versions of the model are available:

- **`ogHANNA`**: The original version (HANNA v1.0.0), trained on binary VLE data (up to 10 bar) and limiting activity coefficients from the Dortmund Data Bank. This version is limited to binary mixtures.
- **`HANNA`** (alias `multHANNA`): The latest version, trained on VLE and LLE data, and applicable to multi-component mixtures.

```@docs
MLThermoProperties.ogHANNA
MLThermoProperties.multHANNA
```

## GRAPPA saturation model

GRAPPA is a graph neural network model for predicting vapor pressures and boiling points of pure components. The model predicts the parameters `A`, `B`, and `C` of the Antoine equation:

```math
\ln(p^s / \text{kPa}) = A - \frac{B}{T / \text{K} + C}
```

On model construction, the Antoine parameters are predicted and a [`SaturationModel`](@ref) is automatically created, which enables the calculation of the vapor pressure via [`saturation_pressure`](@ref) for a given temperature.

```@docs
MLThermoProperties.GRAPPA
```

## ML utilities

The package [`ChemBERTa.jl`](https://github.com/se-schmitt/MLThermoProperties.jl/tree/main/lib/ChemBERTa) contains encoder language models from the ChemBERTa model family.
It is an registered package and can be used independently of `MLThermoProperties.jl`.

```@docs
MLThermoProperties.ChemBERTa.ChemBERTaModel
```
