#to be extended on ClapeyronGlennExt
struct GlennJL{S,T} <: IdealModel
   components::Vector{String}
   species_info::Vector{S}
   intervals::Vector{T}
   reference_state::ReferenceState
   Rgas::Float64
   R0::Vector{Float64}
   references::Vector{String}
end

function __glenn_jl end


"""
    GlennJL <: IdealModel

    GlennJL(calc::Union{Glenn.Calculator,Glenn.ThermoDB}, input;
            Rgas = Clapeyron.Rgas(),
            R0 = 8.314510,
            reference_state = nothing,
            verbose = false,
            strict = true)

## Input arguments

- `calc`: a database of Glenn.jl parameters 
- `input`: input species, it can be strings, integers (ID of the database), `Glenn.SpeciesInfo` or a vector of those elements.
- `Rgas`: molar gas constant used by the model.
- `R0`: molar gas constant used in the fitting of each species. by default it used the original R constant used in the fitting of NASA-7 polynomials.
- `reference_state`: if a change of reference state is needed.
- `verbose`: if set to `true`, displays additional information to the REPL.
- `strict`: if set to `true`, check if all components have the same phase and if their temperature ranges do not intersect


## Description

Ideal model using the NASA-7 polynomial coefficients provided by the `Glenn.jl` package.
The reference state of the NASA-7 polynomials is set such as the enthalpy at p = 1 atm and T = 298.15 K is equal to the formation enthalpy.

## Model Construction Examples
```
# Using the default database
calc = Calculator()

# Using a string
idealmodel = GlennJL(calc,"o2") 

# Using the ID of O2 in the default database
idealmodel = GlennJL(calc,931) 

# Using a SpeciesInfo to search the database
o2 = only(get_available_species(calc, "O2", exact_match = true))
idealmodel = GlennJL(calc,o2) 

# Multiple components, any of the types above can be used instead of strings
idealmodel = GlennJL(calc,["o2","n2"])
```

## Glenn.jl integration

the model `GlennJL` has the following integrations with `Glenn.jl` package:

- `Glenn.calculate_h(model::GlennJL,T,z = [1.0])`
- `Glenn.calculate_s(model::GlennJL,T,z = [1.0])`
- `Glenn.calculate_cp(model::GlennJL,T,z = [1.0])`
- `Glenn.calculate_formation_enthalpy(model::GlennJL,T,z = [1.0])`
- `Glenn.calculate_enthalpy_change(model::GlennJL,T1,T2,z = [1.0])`
- `Glenn.calculate_properties(model::GlennJL,T,z = [1.0])`
- `Glenn.get_properties_range(model::GlennJL,T,z = [1.0])`

Some notes:

- `get_properties` and `get_properties_range` will return NaN for properties outside their ranges. in particular, `get_properties_range` will return a vector of the same size as the input, with `NaN` on invalid inputs.
- dimensionless entropy values returned by `calculate_s` are calculated at p = 1 atm.
```
"""
function GlennJL(a1,a2;kwargs...)
    m = Base.get_extension(Clapeyron,:ClapeyronGlennExt)
    if m === nothing
        error("the `GlennJL` model requires the Package 'Glenn.jl' to be available in the same enviroment as `Clapeyron.jl`. Please add the package if not installed (`]add Glenn`) and then load it in the current enviroment (`using Glenn`).")
    end
    __glenn_jl(a1,a2;kwargs...)
end

paramtype(::Type{GlennJL{S,T}}) where {S,T} = Float64
Rgas(model::GlennJL) = model.Rgas

export GlennJL

