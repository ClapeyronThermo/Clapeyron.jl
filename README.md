![106277775_782656558933856_3832671477006550609_n](/Users/pierrewalker/Downloads/106277775_782656558933856_3832671477006550609_n.jpg)

Welcome to JuliaSAFT! This module intends to provide the variants of the SAFT equation of state, along with the relevant parameters and solvers required to use these equations.

SAFT equations of state currently available:

| EoS         | $A_\mathrm{seg.}$/$A_\mathrm{mono.}$? | $A_\mathrm{chain}$? | $A_\mathrm{assoc.}$? | Parameters? |
| ----------- | ------------------------------------- | ------------------- | -------------------- | ----------- |
| SAFT        |                                       |                     |                      |             |
| PC-SAFT     | :heavy_check_mark:                    | :heavy_check_mark:  |                      |             |
| sPC-SAFT    |                                       |                     |                      |             |
| SAFT-VR SW  |                                       |                     |                      |             |
| soft-SAFT   |                                       |                     |                      |             |
| SAFT-VR Mie | :heavy_check_mark:                    |                     |                      |             |

For group contribution approaches, we provide:

| EoS               | $A_\mathrm{seg.}$/$A_\mathrm{mono.}$? | $A_\mathrm{chain}$? | $A_\mathrm{assoc.}$? | Parameters? |
| ----------------- | ------------------------------------- | ------------------- | -------------------- | ----------- |
| sPC-SAFT          |                                       |                     |                      |             |
| SAFT-$\gamma$ SW  |                                       |                     |                      |             |
| SAFT-$\gamma$ Mie |                                       |                     |                      |             |

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

We will also provide a Tp-flash algorithm (HELD alogrithm).