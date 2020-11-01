[![DOI](https://zenodo.org/badge/267659508.svg)](https://zenodo.org/badge/latestdoi/267659508)

![OpenSAFT_logo](docs/src/assets/logo.svg)

Welcome to OpenSAFT! This module intends to provide the variants of the Statistical Associating Fluid Theory (SAFT) thermodynamic equation of state, along with the relevant parameters and solvers required to use these equations.

Check out the Jupyter notebooks in the ```examples``` directory to see how to set up your model.

SAFT equations of state currently available:

| EoS           | Seg./Mono.?        | Chain?             | Assoc.?            | Parameters?        |
| ------------- | ------------------ | ------------------ | ------------------ | ------------------ |
| SAFT          | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: |
| CK-SAFT       |                    |                    |                    |                    |
| sSAFT         |                    |                    |                    |                    |
| LJ-SAFT       |                    |                    |                    |                    |
| PC-SAFT       | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: |
| sPC-SAFT      | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: |
| SAFT-VR SW    |                    |                    |                    |                    |
| soft-SAFT     |                    |                    |                    |                    |
| SAFT-VR Mie   | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: |
| SAFT-VR Morse |                    |                    |                    |                    |

For group contribution approaches, we provide:

| EoS          | Seg./Mono.?        | Chain? | Assoc.? | Parameters? |
| ------------ | ------------------ | ------ | ------- | ----------- |
| sPC-SAFT     |                    |        |         |             |
| SAFT-*ɣ* SW  |                    |        |         |             |
| SAFT-*ɣ* Mie | :heavy_check_mark: |        |         |             |

We also provide some engineering cubic equations of state for comparison:

| EoS                    | Available?         | Parameters?        |
| ---------------------- | ------------------ | ------------------ |
| van der Waals          | :heavy_check_mark: | :heavy_check_mark: |
| Redlich-Kwong          | :heavy_check_mark: | :heavy_check_mark: |
| Soave-Redlich-Kwong    | :heavy_check_mark: | :heavy_check_mark: |
| Peng-Robinson          | :heavy_check_mark: | :heavy_check_mark: |
| Cubic-Plus-Association |                    |                    |

Properties available:

- Bulk, single-phase properties:

| Property                     | Available?         |
| ---------------------------- | ------------------ |
| Volume                       | :heavy_check_mark: |
| Pressure                     | :heavy_check_mark: |
| Entropy                      | :heavy_check_mark: |
| Internal Energy              | :heavy_check_mark: |
| Enthalpy                     | :heavy_check_mark: |
| Gibbs free energy            | :heavy_check_mark: |
| Helmholtz free energy        | :heavy_check_mark: |
| Isochoric heat capacity      | :heavy_check_mark: |
| Isobaric heat capacity       | :heavy_check_mark: |
| Isentropic compressibility   | :heavy_check_mark: |
| Isothermal compressibility   | :heavy_check_mark: |
| Isobaric (cubic) expansivity | :heavy_check_mark: |
| Speed of sound               | :heavy_check_mark: |
| Joule-Thomson coefficient    | :heavy_check_mark: |

- Two-phase properties:

| Property                  | Available?                         |
| ------------------------- | ---------------------------------- |
| Saturation pressure       | :heavy_check_mark:                 |
| Bubble pressure           | :heavy_check_mark: (SAFT EOS only)​ |
| Dew pressure              |                                    |
| Bubble temperature        |                                    |
| Dew temperature           |                                    |
| Enthalpy of vapourisation | :heavy_check_mark:                 |

- Critical properties (pure components only):

| Property             | Available?         |
| :------------------- | ------------------ |
| Critical temperature | :heavy_check_mark: |
| Critical pressure    | :heavy_check_mark: |
| Critical volume      | :heavy_check_mark: |

We will also provide a Tp-flash algorithm (Rachford-Rice and HELD alogrithm).

Note that at its current stage, OpenSAFT is still in the very early stages of development, and things may be moving around or changing rapidly, but we are very excited to see where this project may go!

# Installing OpenSAFT

OpenSAFT is not yet in the JuliaHub (but it will be soon!).

To load OpenSAFT, launch Julia with

    > julia

Hit the ```]``` key to enter Pkg mode, then type

    Pkg> add git@github.com:ypaul21/OpenSAFT.jl.git
