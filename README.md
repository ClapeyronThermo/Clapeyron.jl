[![DOI](https://zenodo.org/badge/267659508.svg)](https://zenodo.org/badge/latestdoi/267659508)
[![Build Status](https://github.com/ypaul21/Clapeyron.jl/workflows/CI/badge.svg)](https://github.com/ypaul21/Clapeyron.jl/actions)
[![codecov](https://codecov.io/gh/ypaul21/Clapeyron.jl/branch/master/graph/badge.svg?token=ZVGGR4AAFB)](https://codecov.io/gh/ypaul21/Clapeyron.jl)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://ypaul21.github.io/Clapeyron.jl/dev)
[![project chat](https://img.shields.io/badge/zulip-join_chat-brightgreen.svg)](https://julialang.zulipchat.com/#narrow/stream/265161-Clapeyron.2Ejl)

![Clapeyron_logo](docs/Clapeyron_logo.svg)

Welcome to Clapeyron! This module provides both a large library of thermodynamic models and a framework for one to easily implement their own models.

The official manuscript is in ACS. https://pubs.acs.org/doi/10.1021/acs.iecr.2c00326. There is also a preprint available at arxiv: https://arxiv.org/abs/2201.08927

We have also ppresented at the JuliaCon 2021 conference! Feel free to take a look at our talk:

[![Clapeyron.jl: An Extensible Implementation of Equations of State | Paul Yew et al | JuliaCon2021](https://img.youtube.com/vi/Re5qI-9zyIM/0.jpg)](https://www.youtube.com/watch?v=Re5qI-9zyIM "Clapeyron.jl: An Extensible Implementation of Equations of State | Paul Yew et al | JuliaCon2021")

We support many equations of state and properties. Some examples of figures you can create are shown below:

- Isobaric heat capacity of carbon dioxide at 20 MPa:

  ![CO2_cp](docs/CO2_cp.svg) 

- Water VLE envelope:

  ![water_VLE](docs/water_VLE.svg)

- Ethanol+water Pxy diagram at 423.15 K:

  ![ethanol+water](docs/ethanol+water.svg)

- pT-isopleth of methanol+cyclohexane generated using PC-SAFT:

![CH3OH_CyHx](docs/CH3OH_CyHex.svg)

# Installing Clapeyron

To install Clapeyron, launch Julia with

```julia
> julia
```

Hit the ```]``` key to enter Pkg mode, then type

```julia
Pkg> add Clapeyron
```
Or to add the development version:
```julia
Pkg> add https://github.com/ypaul21/Clapeyron.jl#master
```
Exit Pkg mode by hitting backspace.

Now you may begin using functions from the Clapeyron library by entering the command

```
using Clapeyron
```

To remove the package, hit the ```]``` key to enter Pkg mode, then type

```julia
Pkg> rm Clapeyron
```
## Citing Clapeyron

If you are using Clapeyron for your research work, please cite the following:

```
@article{Clapeyron-2022,
    title={Clapeyron.jl: An Extensible, Open-Source Fluid Thermodynamics Toolkit},
    author={Pierre J. Walker, Hon-Wa Yew, and Andrés Riedemann},
    journal={Ind. Eng. Chem. Res.},
    volume={XX},
    number={XX},
    pages={XX--XX},
    year={2022},
    publisher={American Chemical Society},
    doi={doi/10.1021/acs.iecr.2c00326},
    url={https://pubs.acs.org/doi/10.1021/acs.iecr.2c00326}
}
```

## Package in active Development

Note that at its current stage, Clapeyron is still in the early stages of development, and things may be moving around or changing rapidly, but we are very excited to see where this project may go!

We are open to contributions, new models, improved methods and more databases are always appreciated.

If you find any issue, feel free to contact us directly on the Zulip Channel, or open a Github issue. 
