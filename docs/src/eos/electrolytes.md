```@meta
CurrentModule = Clapeyron
```

## Index

```@index
Pages = ["electrolytes.md"]
```

## Explicit vs. Implicit Solvent Models

Electrolyte equations of state model mixtures containing charged species (ions). For such a mixture to be in thermodynamic equilibrium, a fundamental constraint must be imposed: the total electrical charge must sum to zero. This is the condition of **electroneutrality**.

How a model handles charged species defines its approach:

- **Explicit Solvent Models:** These models work directly with individual ions as components (e.g., Na⁺, Cl⁻). While conceptually straightforward at the model level (ions are uniquely defined) performing phase equilibrium calculations requires explicitly satisfying the electroneutrality constraint alongside all other equilibrium conditions (like equality of chemical potentials). This adds a layer of complexity to equilibrium calculations.

- **Implicit Solvent Models:** To circumvent the complexity of the electroneutrality constraint, this approach groups ions into electrically neutral pairs, representing them as a single **salt component** (e.g., NaCl). The remaining components (such as water) are implicitly treated as the solvent. Consequently, any model that operates with salt components rather than individual ions falls under this category.

Algorithms that solve phase equilibria by working directly with explicit solvent models are referred to as using the **Electrochemical Ionic Approach**, whereas algorithms that first transform the problem into an implicit solvent formulation are known to use the **Mean Ionic Approach**.

In Clapeyron.jl, we mainly implement **explicit solvent** models. However, to provide flexibility and connect the two methodologies, we offer the `MeanIonicApproach` model wrapper, which can transform an explicit solvent model into an equivalent implicit solvent representation.

## Main models

```@docs
Clapeyron.ESElectrolyte
Clapeyron.MeanIonicApproach
```

## Ion Models

```@docs
Clapeyron.Born
Clapeyron.DH
Clapeyron.MSA
Clapeyron.DHBorn
Clapeyron.MSABorn
Clapeyron.GCMSABorn
```

## Electrolyte Models

```@docs
Clapeyron.ePCSAFT
Clapeyron.eSAFTVRMie
Clapeyron.SAFTVREMie
Clapeyron.SAFTgammaEMie
```

## Relative Static Permittivity Models

```@docs
Clapeyron.ConstRSP
Clapeyron.LinMixRSP
Clapeyron.Schreckenberg
Clapeyron.ZuoFurst
```