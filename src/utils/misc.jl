using NamedArrays

export create_z

function create_z(model::EoS, z::AbstractArray)
    return NamedArray(z, model.components)
end

"""
    normalize_units(val)
On normal numbers, it is the identity, but on numbers or vectors of `Unitful.Quantity`,it converts the unit to an equivalent SI unit and strips the unit information.

# Examples
```julia-repl
julia> normalize_units(0.0u"°C")
273.15

```
```julia-repl
julia> normalize_units(273.15)
273.15
```
"""
normalize_units(x) =  Unitful.ustrip.(Unitful.upreferred.(x))

