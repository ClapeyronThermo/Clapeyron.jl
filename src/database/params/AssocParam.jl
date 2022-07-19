"""
    AssocParam{T}
Struct holding association parameters.
"""
struct AssocParam{T} <: ClapeyronParam
    name::String
    components::Array{String,1}
    values::Compressed4DMatrix{T,Vector{T}}
    sites::Array{Array{String,1},1}
    sourcecsvs::Array{String,1}
    sources::Array{String,1}
end

function AssocParam(
        name::String,
        components::Vector{String},
        values::MatrixofMatrices,
        allcomponentsites,
        sourcecsvs,
        sources
    ) where T
    _values = Compressed4DMatrix(values)
    return AssocParam(name, components, _values, allcomponentsites, sourcecsvs,sources)
end

function AssocParam(x::AssocParam, name::String = x.name; isdeepcopy = true, sources = x.sources)
    if isdeepcopy
        return AssocParam(
            name,
            x.components,
            deepcopy(x.values),
            x.sites,
            x.sourcecsvs,
            sources
        )
    end
    return AssocParam(
        name,
        x.components,
        x.values,
        x.sites,
        x.sourcecsvs,
        sources
    )
end

function Base.show(io::IO, mime::MIME"text/plain", param::AssocParam{T}) where T
    print(io, "AssocParam{", string(T), "}")
    print(io, param.components)
    println(io, ") with values:")
    comps = param.components
    vals = param.values
    sitenames = param.sites
    for (idx, (i,j), (a,b)) in indices(vals)
        try
        s1 = sitenames[i][a]
        s2 = sitenames[j][b]
        print(io, "(\"", comps[i], "\", \"", s1, "\")")
        print(io, " >=< ")
        print(io, "(\"", comps[j], "\", \"", s2, "\")")
        print(io, ": ")
        println(io, vals.values[idx])
        catch
        println("error at i = $i, j = $j a = $a, b = $b")
        end
    end
end

function Base.show(io::IO, param::AssocParam)
    print(io, typeof(param), "(\"", param.name, "\")")
    print(io, param.values.values)
end

# Operations
function Base.:(+)(param::AssocParam, x::Number)
    values = param.values + x
    return AssocParam(param.name, param.components, values, param.sites ,param.sourcecsvs, param.sources)
end

function Base.:(*)(param::AssocParam, x::Number)
    values = param.values * x
    return AssocParam(param.name, param.components, values, param.sites, param.sourcecsvs, param.sources)
end

function Base.:(^)(param::AssocParam, x::Number)
    values = param.values ^ x
    return AssocParam(param.name, param.components, values, param.sites, param.sourcecsvs, param.sources)
end