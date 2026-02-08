# Installation

This section provides instruction on installing the Julia language and Clapeyron, and provides recommendations on packages users might need to use in conjunction with Clapeyron.

## Installing Julia

The latest version of Julia can be downloaded [here](https://julialang.org/downloads/), with additional instructions specific to the OS used provided [here](https://julialang.org/downloads/platform).
Clapeyron should function on all OS.

If users are unfamiliar with Julia, we recommend some helpful guides to become familiar:

* For novice programmers, there is a short introduction to Julia available [here](https://www.datacamp.com/tutorial/julia-programming-tutorial-for-beginners).
  A longer list of tutorials is available on the Julia [website](https://julialang.org/learning/tutorials/) itself.
* If users are already familiar with MATLAB or Python, there is a great guide available listing the differences between the two [here](https://docs.julialang.org/en/v1/manual/noteworthy-differences/).

For basic usage of Clapeyron, one does not need an in-depth knowledge of Julia.
However, if one wishes to implement their own methods or model, it might be worth familiarising oneself with concepts such as multiple dispatch and broadcasting.

## Installing Clapeyron

Clapeyron.jl is a registered package on JuliaHub.
Installing it can be done with a simple command within the Julia REPL:

```julia
julia> using Pkg

julia> Pkg.add("Clapeyron")
```

If a new version of Clapeyron.jl is released, one can update the package using:

```julia
julia> using Pkg

julia> Pkg.update("Clapeyron")
```

If you want to use the development version, you need to install the `master` branch of the repository, in the following way:

```julia
julia> Pkg.add(url="https://github.com/ClapeyronThermo/Clapeyron.jl", rev="master")
```

## Recommended packages

In order to fully utilise Clapeyron, users may need certain features not included in the package.
Here is a list of packages the developers of Clapeyron recommend using:

* Plotting: The default package for plotting in Julia is [Plots.jl](https://github.com/JuliaPlots/Plots.jl) and can be installed the same way as Clapeyron.
  However, if users are more-familiar with matplotlib, [PyPlot.jl](https://github.com/JuliaPy/PyPlot.jl) is also available but is trickier to install.
* Data storage and manipulation: The default packages in Julia are [DataFrames.jl](https://github.com/JuliaData/DataFrames.jl) and [Tables.jl](https://github.com/JuliaData/Tables.jl).
  Both of these make it easy to store values and then export them into various data types.

## Installing Clapeyron in Python

Clapeyron.jl is also available in python via the [pyclapeyron](github.com/ClapeyronThermo/pyclapeyron) package.
You can install `pyclapeyron` from PyPI via `pip` or `uv`:

```
pip install pyclapeyron
```

or

```
uv add pyclapeyron
```