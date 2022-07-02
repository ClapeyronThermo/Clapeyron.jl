```@meta
CurrentModule = Clapeyron
```
# Clapeyron.jl
A [Julia](http://julialang.org) package for the modelling of fluids using thermodynamic equations of state. These include the standard cubics (van der Waals, Redlich-Kwong, Peng-Robinson, _etc._), SAFT-type equations (PC-SAFT, SAFT-VR Mie, SAFT-$\gamma$ Mie, _etc._), empirical equations (GERG2008, IAWPS95) and many more.

The documentation is laid out as follows:

- **Background**: Find more information about the origin and differences for each equation of state provided in this package, as well as some of the methods used to obtain the various thermodynamic properties.
- **Basic Usage**: Find out how to use the equations of state provided in the package.
- **Customization**: how to use your own parameters, implement your own equation of state or algorithm.
- **Notebook Examples**: A list of available notebooks showcasing different functionalities of the package.
- **Available Models**: A list of all available Equations of State present in the package.
- **API**: A list of all available methods.

### Authors

- [Pierre J. Walker](mailto:pjwalker@caltech.edu), California Institute of Technology
- [Hon-Wa (Paul) Yew](mailto:honwa.yew16@imperial.ac.uk), Imperial College London
- [Andrés Riedemann](mailto:andres.riedemann@gmail.com), University of Concepción

### License

Clapeyron.jl is licensed under the [MIT license](https://github.com/ypaul21/Clapeyron.jl/blob/master/LICENSE.md).

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
    volume={XX},
    number={XX},
    pages={XX--XX},
    year={2022},
    publisher={American Chemical Society},
    doi={doi/10.1021/acs.iecr.2c00326},
    url={https://pubs.acs.org/doi/10.1021/acs.iecr.2c00326}
}
```
## Citing a particular equation of state model

In addition to citing this work, it is encouraged to cite the references to the underlying models used. for that, you can use [`Clapeyron.cite`](@ref) to obtain the references used in a particular model




