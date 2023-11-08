```@meta
CurrentModule = Clapeyron
```
# Clapeyron.jl
An extensible [Julia](http://julialang.org) package for the modelling of fluids using thermodynamic equations of state. These include the standard cubics (van der Waals, Redlich-Kwong, Peng-Robinson, _etc._), SAFT-type equations (PC-SAFT, SAFT-VR Mie, SAFT-$\gamma$ Mie, _etc._), empirical equations (GERG2008, IAWPS95), Activity coefficient models (NRTL, UNIFAC, COSMO-SAC, _etc._) and many more.

The documentation is laid out as follows:

0. **Installation**
1. **Model construction**
1.1 Ideal gas models
1.2 Cubic equations of state
1.3 SAFT equations of state
1.4 Empirical equations of state
1.5 Activity coefficient models
1.6 Composite Models
1.A Available groups
2. **Property Estimation**
2.1 Bulk Properties
2.2 Phase Equilibria
2.3 Flash Calculations
3. **Parameter Estimation**
4. **Extensions**
4.1 CoolProp
4.2 Symbolics
4.3 Measurements
5. **Package Development**
5.1 Developing new methods
5.2 Developing new models
5.3 Contributing to the package
6. **API**

### Authors

- [Pierre J. Walker](mailto:pjwalker@caltech.edu), California Institute of Technology
- [Hon-Wa (Paul) Yew](mailto:honwa.yew16@imperial.ac.uk), Imperial College London
- [Andrés Riedemann](mailto:andres.riedemann@gmail.com), University of Concepción

### License

Clapeyron.jl is licensed under the [MIT license](https://github.com/ClapeyronThermo/Clapeyron.jl/blob/master/LICENSE.md).

### Installation

Clapeyron.jl is a registered package, it can be installed from the general registry by:

```
pkg> add Clapeyron
```

## Citing `Clapeyron.jl`

If you are using Clapeyron for your research work, please cite the following:

```
@article{Clapeyron-2022,
    title={Clapeyron.jl: An Extensible, Open-Source Fluid Thermodynamics Toolkit},
    author={Pierre J. Walker, Hon-Wa Yew, and Andrés Riedemann},
    journal={Ind. Eng. Chem. Res.},
    volume={61},
    number={20},
    pages={7130--7153},
    year={2022},
    publisher={American Chemical Society},
    doi={doi/10.1021/acs.iecr.2c00326},
    url={https://pubs.acs.org/doi/10.1021/acs.iecr.2c00326}
}
```
## Citing a particular equation of state model

In addition to citing this work, it is encouraged to cite the references to the underlying models used. for that, you can use [`Clapeyron.cite`](@ref) to obtain the references used in a particular model




