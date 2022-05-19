abstract type EoSModel end
export EoSModel

groups(model::EoSModel) = model.groups
components(model::EoSModel) = model.components
components(model) = nothing

function eos_length_error(model,l,lz)
    throw(DimensionMismatch("EoS Model has length $l, molar amount vector has length $lz"))
end

function eos_length_check(model,z)
    components(model) === nothing && return nothing  
    lz = length(z)
    lm = length(model)
    if !(lm == lz)
        eos_length_error(model,lm,lz)
    end
    return nothing
end

"""
    eos(model::EoSModel, V, T, z=SA[1.0])

basic Clapeyron function, returns the total Helmholtz free energy.

# Inputs:
- `model::EoSModel` Thermodynamic model to evaluate
- `V` Total volume, in [m³]
- `T` Temperature, in [K]
- `z` mole amounts, in [mol], by default is `@SVector [1.0]`

# Outputs:

- Total Helmholtz free energy, in [J]

by default, it calls `R̄*T*∑(z)*(a_ideal(ideal_model,V,T,z) + a_res(model,V,T,z))` where `ideal_model == idealmodel(model)`, where `a_res` is the reduced residual Helmholtz energy and `a_ideal` is the reduced ideal Helmholtz energy.
You can mix and match ideal models if you provide:
- `idealmodel(model)`: extracts the ideal model from your Thermodynamic model
- `a_res(model,V,T,z)`: residual reduced Helmholtz free energy

"""
function eos(model::EoSModel, V, T, z=SA[1.0])
    eos_length_check(model,z)
    return N_A*k_B*sum(z)*T * (a_ideal(idealmodel(model),V,T,z)+a_res(model,V,T,z))
end
"""
    idealmodel(model::EoSModel)
    
retrieves the ideal model from the input's model.

# Examples:

```julia-repl
julia> pr = PR(["water"],idealmodel=IAPWS95Ideal)   
PR{IAPWS95Ideal} with 1 component:
 "water"
Contains parameters: a, b, acentricfactor, Tc, Mw   
julia> Clapeyron.idealmodel(pr)
IAPWS95Ideal()
```

"""
idealmodel(model::EoSModel) = model.idealmodel

"""
    eos_res(model::EoSModel, V, T, z=SA[1.0])

basic Clapeyron function, returns the residual Helmholtz free energy.

# Inputs:
- `model::EoSModel` Thermodynamic model to evaluate
- `V` Total volume, in [m³]
- `T` Temperature, in [K]
- `z` mole amounts, in [mol], by default is `@SVector [1.0]`

# Outputs:

- Residual Helmholtz free energy, in [J]

by default, it calls `R̄*T*∑(z)*(a_res(model,V,T,z))` where `a_res` is the reduced residual Helmholtz energy.
"""
function eos_res(model::EoSModel, V, T, z=SA[1.0])
    eos_length_check(model,z)
    return N_A*k_B*sum(z)*T*(a_res(model,V,T,z))
end

Base.broadcastable(model::EoSModel) = Ref(model)

"""
    @comps

This macro is an alias to

    1:length(model)

The caveat is that `model` has to exist in the local namespace.
"""
macro comps()
    return quote
        1:length(model)
    end |> esc
end

has_sites(::Type{<:EoSModel}) = false
has_groups(::Type{<:T}) where T = false
has_sites(::T) where T = has_sites(T)
has_groups(::T) where T = has_groups(T)

