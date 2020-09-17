![JuliaSAFT_logo](./docs/JuliaSAFT_logo.jpg)

Welcome to JuliaSAFT! This module intends to provide the variants of the SAFT equation of state, along with the relevant parameters and solvers required to use these equations.

SAFT equations of state currently available:

| EoS           | Seg./Mono.?        | Chain?             | Assoc.? | Parameters? |
| ------------- | ------------------ | ------------------ | ------- | ----------- |
| SAFT          | :heavy_check_mark: | :heavy_check_mark: |         |             |
| CK-SAFT       |                    |                    |         |             |
| sSAFT         |                    |                    |         |             |
| LJ-SAFT       |                    |                    |         |             |
| PC-SAFT       | :heavy_check_mark: | :heavy_check_mark: |         |             |
| sPC-SAFT      | :heavy_check_mark: | :heavy_check_mark: |         |             |
| SAFT-VR SW    |                    |                    |         |             |
| soft-SAFT     |                    |                    |         |             |
| SAFT-VR Mie   | :heavy_check_mark: | :heavy_check_mark: |         |             |
| SAFT-VR Morse |                    |                    |         |             |

For group contribution approaches, we provide:

| EoS               | Seg./Mono.? | Chain? | Assoc.? | Parameters? |
| ----------------- | ----------- | ------ | ------- | ----------- |
| sPC-SAFT          |             |        |         |             |
| SAFT-$\gamma$ SW  |             |        |         |             |
| SAFT-$\gamma$ Mie |             |        |         |             |

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

| Property                  | Available?         |
| ------------------------- | ------------------ |
| Saturation pressure       | :heavy_check_mark: |
| Bubble pressure           |                    |
| Dew pressure              |                    |
| Bubble temperature        |                    |
| Dew temperature           |                    |
| Enthalpy of vapourisation | :heavy_check_mark: |

- Critical properties (pure components only):

| Property             | Available?         |
| :------------------- | ------------------ |
| Critical temperature | :heavy_check_mark: |
| Critical pressure    | :heavy_check_mark: |
| Critical volume      | :heavy_check_mark: |

We will also provide a Tp-flash algorithm (HELD alogrithm).