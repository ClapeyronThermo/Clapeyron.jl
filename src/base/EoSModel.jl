abstract type EoSModel end
export EoSModel

"""
    eos(model::EoSModel, V, T, z=SA[1.0])

Returns the total Helmholtz free energy.

# Inputs:
- `model::EoSModel` Thermodynamic model to evaluate
- `V` Total volume, in [m³]
- `T` Temperature, in [K]
- `z` mole amounts, in [mol], by default is `@SVector [1.0]`

# Outputs:
- Total Helmholtz free energy, in [J]

by default, it calls `R̄*T*∑(z)*(a_ideal(ideal_model,V,T,z) + a_res(model,V,T,z))` where `ideal_model == idealmodel(model)`, where `a_res` is the reduced residual Helmholtz energy and `a_ideal` is the reduced ideal Helmholtz energy.
You can mix and match ideal models if you provide:
- `[idealmodel](@ref)(model)`: extracts the ideal model from your Thermodynamic model
- `[a_res](@ref)(model,V,T,z)`: residual reduced Helmholtz free energy

"""
function eos(model::EoSModel, V, T, z=SA[1.0])
    return N_A*k_B*sum(z)*T * (a_ideal(idealmodel(model),V,T,z)+a_res(model,V,T,z))
end
"""
    idealmodel(model::EoSModel)
    
retrieves the ideal model from the input's model.

# Examples:

```julia-repl
julia> pr = PR(["water"],idealmodel=MonomerIdeal)
PR{MonomerIdeal, PRAlpha, NoTranslation, vdW1fRule} with 1 component:
 "water"
Contains parameters: a, b, Tc, Pc, Mw 
julia> Clapeyron.idealmodel(pr)
MonomerIdeal with 1 component:
 "water"
Contains parameters: Mw
```
"""
idealmodel(model::EoSModel) = model.idealmodel

"""
    eos_res(model::EoSModel, V, T, z=SA[1.0])

Returns the residual Helmholtz free energy.

# Inputs:
- `model::EoSModel` Thermodynamic model to evaluate
- `V` Total volume, in [m³]
- `T` Temperature, in [K]
- `z` mole amounts, in [mol], by default is `@SVector [1.0]`

# Outputs:
- Residual Helmholtz free energy, in [J]

by default, it calls `R̄*T*∑(z)*(a_res(model,V,T,z))` where [`a_res`](@ref) is the reduced residual Helmholtz energy.
"""
function eos_res(model::EoSModel, V, T, z=SA[1.0])
    return N_A*k_B*sum(z)*T*(a_res(model,V,T,z))
end


"""
    a_res(model::EoSModel, V, T, z,args...)

Reduced residual Helmholtz free energy.

# Inputs:
- `model::EoSModel` Thermodynamic model to evaluate
- `V` Total volume, in [m³]
- `T` Temperature, in [K]
- `z` mole amounts, in [mol], by default is `@SVector [1.0]`

# Outputs:
- Residual Helmholtz free energy, no units

You can define your own EoS by adding a method to `a_res` that accepts your custom model. 
"""
function a_res end
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

This list will displayed by each `EoSModel` on future versions. you can enable/disable this by setting `ENV["CLAPEYRON_SHOW_REFERENCES"] = "TRUE"/"FALSE"`
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
    