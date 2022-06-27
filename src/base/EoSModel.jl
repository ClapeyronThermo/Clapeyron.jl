abstract type EoSModel end
export EoSModel

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
    return N_A*k_B*sum(z)*T*(a_res(model,V,T,z))
end

Base.broadcastable(model::EoSModel) = Ref(model)

"""
    @comps

This macro is an alias to

    1:length(model)

The caveat is that `model` has to exist in the local namespace.
`model` is expected to be an EoSModel type that contains the `icomponents` field.
`icomponents` is an iterator that goes through all component indices.
"""
macro comps()
    return quote
        1:length(model)
    end |> esc
end

has_sites(::Type{<:EoSModel}) = false
has_groups(::Type{<:EoSModel}) = false
has_sites(::T) where T<:EoSModel = has_sites(T)
has_groups(::T) where T<:EoSModel = has_groups(T)

"""
    doi(model)

Returns a Vector of strings containing the top-level bibliographic references of the model, in DOI format.

```julia-repl
julia> umr = UMRPR(["water"],idealmodel = WalkerIdeal);Clapeyron.doi(umr)
1-element Vector{String}:
 "10.1021/I160057A011"
```
"""
function doi(model)
    if hasfield(typeof(model),:references)
        return model.references
    else
        return String[]
    end
end

"""
    cite(model)

Returns a Vector of strings containing all bibliographic references of the model, in DOI format. this includes any nested models.

```julia-repl
julia> umr = UMRPR(["water"],idealmodel = WalkerIdeal);Clapeyron.cite(umr) #should cite UMRPR, UNIFAC, WalkerIdeal
4-element Vector{String}:
 "10.1021/I160057A011"
 "10.1021/ie049580p"
 "10.1021/i260064a004"
 "10.1021/acs.jced.0c00723"
```

Starting from Clapeyron 0.3.7, this list is displayed by each `EoSModel`. you can disable this by setting `ENV["CLAPEYRON_SHOW_REFERENCES"] = "FALSE"`
"""
function cite(model::EoSModel)
    keys = fieldnames(typeof(model))
    res = doi(model)
    
    for key in keys
        val = getfield(model,key)
        if val isa EoSModel
            append!(res,cite(val))
        elseif key == :references
            append!(res,val)
        end
    end
    return unique!(res)
end
    