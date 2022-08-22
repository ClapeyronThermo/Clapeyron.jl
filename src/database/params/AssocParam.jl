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

function Base.copyto!(dest::AssocParam,src::AssocParam) #used to set params
    #key check
    dest.components == src.components || throw(DimensionMismatch("components of source and destination single parameters are not the same for $dest"))
    copyto!(dest.values,src.values)
    dest_sites = dest.sites
    src_sites = src.sites
    for (i,site) in enumerate(dest_sites)
        copy!(site,src_sites[i]) #copy also changes size
    end
    return dest
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
    l = length(param.values.values)
    print(io, ") with ", l, " value",ifelse(l==1,"","s"),":")
    l != 0 && println(io)
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
        print(io, vals.values[idx])
        l != idx && println(io)
        catch
        println("error at i = $i, j = $j a = $a, b = $b")
        end
    end
end

function Base.show(io::IO, param::AssocParam)
    print(io, typeof(param), "(\"", param.name, "\")")
    print(io, param.values.values)
end