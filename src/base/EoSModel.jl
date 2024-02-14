abstract type EoSModel end

function a_eos(model::EoSModel, V, T, z=SA[1.0])
    maybe_ideal = idealmodel(model)
    ideal = maybe_ideal !== nothing ? maybe_ideal : model
    return a_ideal(ideal,V,T,z) + a_res(model,V,T,z) + reference_state_eval(model,V,T,z)
end

"""
    Rgas(model)
    Rgas()

Returns the gas constant used by an `EoSModel`.

By default, uses the current 2019 definition: `R̄` = 8.31446261815324 [J⋅K⁻¹⋅mol⁻¹]. you can call `Rgas()` to obtain this value.

"""
Rgas(model) = R̄
Rgas() = R̄

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
    return Rgas(model)*sum(z)*T * a_eos(model,V,T,z)
end

"""
    idealmodel(model::EoSModel)

retrieves the ideal model from the input's model. if the model is already an idealmodel, return `nothing`
# Examples:
```julia-repl
julia> pr = PR(["water"],idealmodel = MonomerIdeal)
PR{MonomerIdeal, PRAlpha, NoTranslation, vdW1fRule} with 1 component:
 "water"
Contains parameters: a, b, Tc, Pc, Mw
julia> ideal = idealmodel(pr)
MonomerIdeal with 1 component:
 "water"
Contains parameters: Mw
julia> idealmodel(ideal) == nothing
true
```
"""
idealmodel(model::EoSModel) = __idealmodel(model::EoSModel)

@generated function __idealmodel(model::EoSModel)
    if hasfield(model,:idealmodel)
        return :(getfield(model,:idealmodel))
    else
        return :(nothing)
    end
end

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
    return Rgas(model)*sum(z)*T*a_res(model,V,T,z)
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
Base.transpose(model::EoSModel) = model
Base.eltype(model::EoSModel) = __eltype(model)
@pure function __eltype(model::T) where T <:EoSModel
    if hasfield(T,:params)
        return eltype(model.params)
    else
        return Float64
    end
end

Base.summary(model::EoSModel) = string(parameterless_type(model))

"""
    @comps

This macro is an alias to `1:length(model)`
The caveat is that `model` has to exist in the local namespace.
`model` is expected to any struct that has length defined in terms of the amount of components.
"""
macro comps()
    return quote
        1:length(model)
    end |> esc
end

"""
    R̄

This macro is an alias to `Rgas(model)`
The caveat is that `model` has to exist in the local namespace.
"""
macro R̄()
    return quote
        Rgas(model)
    end |> esc
end

has_sites(::T) where T <: EoSModel = has_sites(T)
has_sites(::Type{T}) where T <: EoSModel = _has_sites(T)

@pure function _has_sites(::Type{T}) where T <: EoSModel
    s1 = hasfield(T,:sites)
    if s1
       return fieldtype(T,:sites) == SiteParam
    end
    return false
end

has_groups(::T) where T <: EoSModel = has_groups(T)
has_groups(::Type{T}) where T <: EoSModel = _has_groups(T)

@pure function _has_groups(::Type{T}) where T <: EoSModel
    s1 = hasfield(T,:groups)
    if s1
       return fieldtype(T,:groups) <: GroupParameter
    end
    return false
end

Base.length(model::T) where T <:EoSModel = _eos_length(model,Val(has_groups(model)))

_eos_length(model::EoSModel,::Val{true}) = length(model.groups.components)
_eos_length(model::EoSModel,::Val{false}) = length(model.components)

"""
    doi(model)
Returns a Vector of strings containing the top-level bibliographic references of the model, in DOI format. if there isn't a `references` field, it defaults to `default_references(model)`
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
        return default_references(model)
    end
end
"""
    default_references(::Type{<:EoSModel})::Vector{String}

Return the default references of a model. If you are using the `@newmodel`, `@newmodelsimple` or `@newmodelgc` macros, define this function to set the references for the defined EoS.

"""
function default_references end
default_references(m::EoSModel) = default_references(parameterless_type(m))
default_references(M) = String[]

"""
    cite(model,out = :doi)

Returns a Vector of strings containing all bibliographic references of the model, in the format indicated by the `out` argument. this includes any nested models.

```julia-repl
julia> umr = UMRPR(["water"],idealmodel = WalkerIdeal);Clapeyron.cite(umr) #should cite UMRPR, UNIFAC, WalkerIdeal
4-element Vector{String}:
 "10.1021/I160057A011"
 "10.1021/ie049580p"
 "10.1021/i260064a004"
 "10.1021/acs.jced.0c00723"
```
the `out` argument supports two values:
- `:doi`: returns the stored values on each EoS. by default those are DOI identifiers.
- `:bib`: returns BibTeX entries. to use this, an internet connection is required.

```julia-repl
julia> model = SAFTVRQMie(["helium"])
SAFTVRQMie{BasicIdeal} with 1 component:
 "helium"
Contains parameters: Mw, segment, sigma, lambda_a, lambda_r, epsilon

julia> Clapeyron.cite(model,:bib)
2-element Vector{String}:
 "@article{Aasen_2019,\n\tdoi = {10" ⋯ 463 bytes ⋯ "Journal of Chemical Physics}\n}"
 "@article{Aasen_2020,\n\tdoi = {10" ⋯ 452 bytes ⋯ "Journal of Chemical Physics}\n}"
```

This list will displayed by each `EoSModel` on future versions. you can enable/disable this by setting `ENV["CLAPEYRON_SHOW_REFERENCES"] = "TRUE"/"FALSE"`
"""
function cite(model::EoSModel,out = :doi)
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
    unique!(res)
    if out == :doi
        return res
    elseif out == :bib
        return doi2bib.(res)
    else
        error("invalid out value $(out)")
    end
end

"""
    recombine!(model::EoSModel)
Recalculate all mixing rules, combining rules and parameter caches inside an `EoSModel`.
"""
function recombine! end

function setreferences!(model,references)
    oldrefs = model.references
    resize!(oldrefs,length(references))
    oldrefs .= references
end

export EoSModel, eos, has_groups, has_sites, Rgas
