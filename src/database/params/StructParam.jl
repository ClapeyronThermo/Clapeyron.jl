"""
    StructParam{T}
Struct holding structure parameters.
"""
struct StructParam{T} <: ClapeyronParam
    name::String
    components::Array{String,1}
    values::Compressed3DMatrix{T,Vector{T}}
    sites::Array{Array{String,1},1}
    sourcecsvs::Array{String,1}
    sources::Array{String,1}
end

function StructParam(
        name::String,
        components::Vector{String},
        values::MatrixofMatrices,
        allcomponentsites,
        sourcecsvs,
        sources
    ) where T
    _values = Compressed3DMatrix(values)
    return StructParam(name, components, _values, allcomponentsites, sourcecsvs,sources)
end

function StructParam(x::StructParam, name::String = x.name; isdeepcopy = true, sources = x.sources)
    if isdeepcopy
        return StructParam(
            name,
            x.components,
            deepcopy(x.values),
            x.sites,
            x.sourcecsvs,
            sources
        )
    end
    return StructParam(
        name,
        x.components,
        x.values,
        x.sites,
        x.sourcecsvs,
        sources
    )
end

function Base.show(io::IO, mime::MIME"text/plain", param::StructParam{T}) where T
    print(io, "StructParam{", string(T), "}")
    print(io, param.components)
    println(io, ") with values:")
    comps = param.components
    vals = param.values
    sitenames = param.sites
    for (idx, i, (k,l)) in indices(vals)
        try
        g1 = sitenames[i][k]
        g2 = sitenames[i][l]
        print(io, "(\"", comps[i], "\", \"", g1, "\")")
        print(io, " - ")
        print(io, "(\"", g2, "\")")
        print(io, ": ")
        println(io, vals.values[idx])
        catch
        println("error at i = $i, k = $k, l = $l")
        end
    end
end

function Base.show(io::IO, param::StructParam)
    print(io, typeof(param), "(\"", param.name, "\")")
    print(io, param.values.values)
end

# Operations
function Base.:(+)(param::StructParam, x::Number)
    values = param.values + x
    return StructParam(param.name, param.components, values, param.sites ,param.sourcecsvs, param.sources)
end

function Base.:(*)(param::StructParam, x::Number)
    values = param.values * x
    return StructParam(param.name, param.components, values, param.sites, param.sourcecsvs, param.sources)
end

function Base.:(^)(param::StructParam, x::Number)
    values = param.values ^ x
    return StructParam(param.name, param.components, values, param.sites, param.sourcecsvs, param.sources)
end