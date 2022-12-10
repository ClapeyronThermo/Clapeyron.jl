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

Base.eltype(param::AssocParam{T}) where T = T

function Base.getindex(param::AssocParam,i::Int) 
    Base.checkbounds(param.components,i)
    getindex(param.values,i,i)
end

function Base.getindex(param::AssocParam,i::Int,j::Int) 
    Base.checkbounds(param.components,max(i,j))
    getindex(param.values,i,j)
end


function AssocParam(
        name::String,
        components::Vector{String},
        values::MatrixofMatrices,
        allcomponentsites,
        sourcecsvs,
        sources
    )
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

#convert utilities
function Base.convert(::Type{AssocParam{Float64}},param::AssocParam{Int})
    assoc_values = param.values
    new_assoc_values = Float64.(assoc_values.values)
    values = Compressed4DMatrix(new_assoc_values,assoc_values.outer_indices,assoc_values.inner_indices,assoc_values.outer_size,assoc_values.inner_size)
    return AssocParam(param.name,param.components,values,param.sites,param.sourcecsvs,param.sources)
end

function Base.convert(::Type{AssocParam{Bool}},param::AssocParam{<:Union{Int,Float64}})
    assoc_values = param.values
    @assert all(z->(isone(z) | iszero(z)),assoc_values.values)
    new_assoc_values = Array(Bool.(assoc_values.values))
    values = Compressed4DMatrix(new_assoc_values,assoc_values.outer_indices,assoc_values.inner_indices,assoc_values.outer_size,assoc_values.inner_size)

    return AssocParam(param.name,param.components,values,param.sites,param.sourcecsvs,param.sources)
end

function Base.convert(::Type{AssocParam{Int}},param::AssocParam{Float64})
    assoc_values = param.values
    @assert all(z->isinteger(z),assoc_values.values)
    new_assoc_values = Int.(assoc_values.values)
    values = Compressed4DMatrix(new_assoc_values,assoc_values.outer_indices,assoc_values.inner_indices,assoc_values.outer_size,assoc_values.inner_size)

    return AssocParam(param.name,param.components,values,param.sites,param.sourcecsvs,param.sources)
end

function Base.convert(::Type{AssocParam{String}},param::AssocParam{<:AbstractString})
    assoc_values = param.values
    new_assoc_values = String.(assoc_values.values)
    values = Compressed4DMatrix(new_assoc_values,assoc_values.outer_indices,assoc_values.inner_indices,assoc_values.outer_size,assoc_values.inner_size)
    return AssocParameter(param.name,param.components,values,param.sites,param.sourcecsvs,param.sources)
end

#trying to break stack overflow on julia 1.6
function Base.convert(::Type{AssocParam{String}},param::AssocParam{String})
    return param
end