# Clapeyron.jl
A [Julia](http://julialang.org) package for the modelling of fluids using thermodynamic equations of state. These include the standard cubics (van der Waals, Redlich-Kwong, Peng-Robinson, _etc._), SAFT-type equations (PC-SAFT, SAFT-VR Mie, SAFT-$\gamma$ Mie, _etc._), empirical equations (GERG2008, IAWPS95) and many more.

The documentation is laid out as follows:

- Background: Find more information about the origin and differences for each equation of state provided in this package, as well as some of the methods used to obtain the various thermodynamic properties
- User guide: Find out how to use the equations of state provided in the package, how to use your own parameters, implement your own equation of state or algorithm.

### Authors

- [Pierre J. Walker](mailto:pjwalker@caltech.edu), California Institute of Technology
- [Hon-Wa (Paul) Yew](mailto:honwa.yew16@imperial.ac.uk), Imperial College London
- [Andrés Riedemann](mailto:andres.riedemann@gmail.com), University of Concepción

### License

Clapeyron.jl is licensed under the [MIT license](https://github.com/ypaul21/Clapeyron.jl/blob/master/LICENSE.md).

### Installation

Clapeyron.jl is a registered package, it can be installed from the general registry by:

```julia
pkg> add Clapeyron
```



