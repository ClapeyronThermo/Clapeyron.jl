const FIJ_TYPE = Clapeyron.PairParameter{Float64, SparseArrays.SparseMatrixCSC{Float64, Int64}}

#= in structs.jl
struct EmpiricDepartureValues <: MultiParameterParam
    polexpgauss::PolExpGaussTerm
    F::Float64
end =#

__type_string(ℙ::Type{EmpiricDepartureValues}) = "Departure"

#for showing in SparseMatrix context.
Base.zero(x::EmpiricDepartureValues) = zero(typeof(x))
Base.zero(::Type{EmpiricDepartureValues}) = EmpiricDepartureValues(PolExpGaussTerm(),0.0)
Base.iszero(x::EmpiricDepartureValues) = iszero(x.F)

#for deleting a departure model
Base.convert(::Type{EmpiricDepartureValues},::Nothing) = zero(EmpiricDepartureValues)
function Base.show(io::IO,x::EmpiricDepartureValues) 
    print(io,"aᵢⱼ(")
    Fij = x.F
    Fij != 0 && print(io,"F = $Fij, ")
    term = x.polexpgauss
    k_pol,k_exp,k_gauss = term.iterators
    l_pol,l_exp,l_gauss = length(k_pol),length(k_exp),length(k_gauss)
    text = String[]
    vals = Int[]
    
    l_pol != 0 && (push!(text,"pow"),push!(vals,l_pol))
    l_exp != 0 &&  (push!(text,"exp"),push!(vals,l_exp))
    l_gauss != 0 && (push!(text,"gauss"),push!(vals,l_gauss))
    show_pairs(io,text,vals," = ",quote_string = false,pair_separator = ", ",)
    print(io,")")
end

"""
    create_departure(data,F = 1.0;verbose = false)

Creates a departure model for use in a `MultiFluid` model with `EmpiricDeparture`.

If `data` is a `String` and starts with `{` or `[`, it will be recognized as JSON text. The text will be parsed as a file location otherwise. 
You can pass a `Dict` or `NamedTuple` if you want to skip the JSON parsing.

## Examples

```julia
d1 = create_departure("/data/EthanePropane.json",0.9) # reading from a file

dep = Dict(
    #reduced CoolProp Departure format, you only need the type and parameters.
    #ResidualHelmholtzPower would work too.
    :type => "Exponential",
    :n => [1,1,1,1],
    :t => [1,1,1,1],
    :d => [1,1,1,1],
    :l => [1,1,1,1],
)

d2 = create_departure(dep) #F is set to 1.0
```
"""
function create_departure(x,F = 1.0;verbose = false)
    X = typeof(x)
    if X <: AbstractString
        if first(x) in ('{','[') #inline JSON
            data = JSON.parse(x; dicttype=Dict{Symbol,Any})
            return create_departure(data,F;verbose = verbose)
        else #location
            paths = flattenfilepaths(String[],x)
            length(paths) != 1 && throw(error("multiple files detected for $x"))
            json_string = read(only(path), String)
            data = JSON.parse(json_string; dicttype=Dict{Symbol,Any})
            return create_departure(data,F;verbose = verbose)
        end 
    else
        return _parse_residual(EmpiricDepartureValues,x;verbose = verbose,Fij = F)
    end
    
end

struct EmpiricDepartureParam <: EoSParam
    F::PairParam{Float64}
    parameters::PairParameter{EmpiricDepartureValues, SparseArrays.SparseMatrixCSC{EmpiricDepartureValues, Int64}}
end

@newmodelsimple EmpiricDeparture MultiFluidDepartureModel EmpiricDepartureParam

"""
EmpiricDeparture <: MultiFluidDepartureModel
    EmpiricDeparture(components;
    userlocations = String[],
    verbose = false)

## Input parameters
none
- `F`: Pair Parameter (`Float64`) - Binary Interaction Parameter (no units).
- `parameters`: Pair Parameter (`String`) - JSON data containing the departure terms for the binary pair.

## Description

Departure that uses empiric departure functions:

```
aᵣ = ∑xᵢaᵣᵢ(δ,τ) + Δa
Δa = ∑xᵢxⱼFᵢⱼaᵣᵢⱼ(δ,τ)

aᵣᵢⱼ = ∑nᵢⱼ₋ₖδ^(dᵢⱼ₋ₖ)*τ^(tᵢⱼ₋ₖ) +
    ∑nᵢⱼ₋ₖδ^(dᵢⱼ₋ₖ)τ^(tᵢⱼ₋ₖ)*exp(-gᵢⱼ₋ₖδ^lᵢⱼ₋ₖ) +
    ∑nᵢⱼ₋ₖδ^(dᵢⱼ₋ₖ)τ^(tᵢⱼ₋ₖ)*exp(ηᵢⱼ₋ₖ(δ-εᵢⱼ₋ₖ)^2 + βᵢⱼ₋ₖ(τ-γᵢⱼ₋ₖ)^2)

```

"""
EmpiricDeparture
default_locations(::Type{EmpiricDeparture}) = ["Empiric/departure/empiric_unlike.csv"]
default_getparams_arguments(::Type{EmpiricDeparture},userlocations,verbose) = ParamOptions(;userlocations,verbose,asymmetricparams = ["F","parameters"])
function transform_params(::Type{EmpiricDeparture},params,components,verbose)
    raw_parameters = params["parameters"]
    F = params["F"]
    s1,s2 = size(F.values)
    𝕊 = sparse((!).(raw_parameters.ismissingvalues))
    dropzeros!(𝕊)
    parsed_parameters = SparseMatrixCSC{EmpiricDepartureValues,Int}(𝕊.m, 𝕊.n, 𝕊.colptr, 𝕊.rowval, similar(𝕊.nzval,EmpiricDepartureValues))
    #parse JSON string to create EmpiricDepartureValues
    for i in 1:s1
        for j in 1:s2
            if !raw_parameters.ismissingvalues[i,j]
                Fij = F[i,j]
                if !iszero(Fij)
                    parsed_parameters[i,j] = create_departure(raw_parameters[i,j],Fij;verbose = verbose)
                else
                    #raw_parameters.ismissingvalues[i,j] = true
                end
            end
        end
    end
    #compress
    parameters = PairParameter(raw_parameters.name,components,parsed_parameters,raw_parameters.ismissingvalues,raw_parameters.sourcecsvs,raw_parameters.sources)
    params["parameters"] = parameters
    return params
end

function multiparameter_a_res(model::MultiFluid,V,T,z,departure::EmpiricDeparture,δ,τ,∑z = sum(z))
    lnδ = log(δ)
    lnτ = log(τ)
    aᵣ = multiparameter_a_res0(model,V,T,z,δ,τ,lnδ,lnτ,∑z)
    _0 = zero(aᵣ)
    ℙ = departure.params.parameters.values
    ℙ_nonzeros = nonzeros(ℙ)
    rows = rowvals(ℙ)
    isone(length(z)) && return aᵣ
    iszero(nnz(ℙ)) && return aᵣ
    Δa = zero(aᵣ)
    @inbounds for j ∈ @comps
        zⱼ = z[j]
        for ii ∈ nzrange(ℙ, j)
            i = rows[ii]
            ℙᵢⱼ = ℙ_nonzeros[ii]
            Fᵢⱼ = ℙᵢⱼ.F
            iszero(Fᵢⱼ) && continue
            aᵢⱼ = reduced_a_res(ℙᵢⱼ,δ,τ,lnδ,lnτ)
            Δa +=z[i]*zⱼ*Fᵢⱼ*aᵢⱼ
        end
     end
    return aᵣ + Δa/(∑z*∑z)
end

function reduced_a_res(ℙ::EmpiricDepartureValues,δ,τ,lnδ = log(δ),lnτ = log(τ))
    return a_term(ℙ.polexpgauss,δ,τ,lnδ,lnτ,zero(δ+τ))
end


"""
    departure_functions(model::MultiFluid)

If the model is using a `EmpiricDeparture` departure model, return the matrix of departure functions. You can set a departure in the following way:

```
using CoolProp #load CoolProp models
model = MultiFluid(["helium","methanol"],mixing = LorentzBerthelotMixing)
dep_mat = departure_functions(model)

dep = Dict(
    #reduced CoolProp Departure format, you only need the type and parameters.
    #ResidualHelmholtzPower would work too.
    type => "Exponential",
    n => [1,1,1,1],
    t => [1,1,1,1],
    d => [1,1,1,1],
    l => [1,1,1,1],
)


dep_mat[1,2] = create_departure(dep,F)

#if you want to delete a departure model:

dep_mat[1,2] = nothing
using SparseArrays
```

"""
function departure_functions(model::MultiFluid)
    dep = model.departure
    dep isa EmpiricDeparture || throw(error("invalid departure model. espected `EmpiricDeparture`, got $(typeof(dep))"))
    return dep.params.parameters.values
end

export EmpiricDeparture, departure_functions, create_departure
