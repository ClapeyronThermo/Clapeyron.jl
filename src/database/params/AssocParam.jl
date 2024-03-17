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

function AssocParam(name::String,components::Vector{String},vals::Compressed4DMatrix{T},sites = nothing,sourcecsvs = String[], sources = String[]) where T
    _sites = if sites === nothing
        ss = [String[] for _ in 1:length(components)]
        for (idx, (i,j), (a,b)) in indices(vals)
            s_i = ss[i]
            s_j = ss[j]
            if length(s_i) < a
                resize!(s_i,a)
            end

            if length(s_j) < b
                resize!(s_j,b)
            end
            s_i[a] = "$(components[i])/site $a"
            s_j[b] = "$(components[j])/site $b"
        end
        ss
    else
        sites
    end

    vals_length = maximum(maximum,vals.outer_indices)
    param_length_check(AssocParam,name,length(components),vals_length)

    return AssocParam{T}(name,components,vals,_sites,sourcecsvs,sources)
end

function AssocParam(name::String,components::Vector{String})
    sites = [String[] for i in 1:length(components)]
    values = Compressed4DMatrix{Float64}()
    AssocParam{Float64}(name,components,values,sites,String[],String[])
end

function Base.copyto!(dest::AssocParam,src::Base.Broadcast.Broadcasted)
    Base.copyto!(dest.values.values,src)
    return dest
end

Base.broadcastable(param::AssocParam) = param.values.values
Base.BroadcastStyle(::Type{<:AssocParam}) = Broadcast.Style{AssocParam}()


function Base.copyto!(dest::AssocParam,src::AbstractArray)
    Base.copyto!(dest.values.values,src)
    return dest
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

Base.size(param::AssocParam) = size(param.values.values)

function Base.getindex(param::AssocParam,i::Int) 
    Base.checkbounds(param.components,i)
    getindex(param.values,i,i)
end

function Base.getindex(param::AssocParam,i::AbstractString)
    idx = findfirst(isequal(i),param.components)
    isnothing(idx) && throw(BoundsError(param,-1))
    getindex(param.values,idx,idx)
end

function Base.getindex(param::AssocParam,i::Int,j::Int) 
    Base.checkbounds(param.components,max(i,j))
    getindex(param.values,i,j)
end

function Base.getindex(param::AssocParam,i::AbstractString,j::AbstractString)
    idx_i = findfirst(isequal(i),param.components)
    idx_j = findfirst(isequal(j),param.components)
    if idx_i === nothing && idx_j === nothing
        throw(BoundsError(param,(-1,-1)))
    elseif idx_i === nothing
        throw(BoundsError(param,(-1,idx_j)))
    elseif idx_j === nothing
        throw(BoundsError(param,(idx_i,-1)))
    else
       return param[idx_i::Int,idx_j::Int]
    end
end

function Base.getindex(param::AssocParam,i::NTuple{2,String},j::NTuple{2,String})
    ii = first(i)
    jj = first(j)
    aa = last(i)
    bb = last(j)
    idx_i,idx_j = _str_to_idx(param,ii,jj)
    sites_i,sites_j = param.sites[idx_i],param.sites[idx_j]
    _idx_a = findfirst(isequal(aa),sites_i)
    _idx_b = findfirst(isequal(bb),sites_j)
    idx_a = _idx_a === nothing ? 0 : _idx_a
    idx_b = _idx_b === nothing ? 0 : _idx_b
    param[idx_i::Int,idx_j::Int][idx_a::Int,idx_b::Int]
end

function Base.setindex!(param::AssocParam,val,i::NTuple{2,String},j::NTuple{2,String})
    ii = first(i)
    jj = first(j)
    aa = last(i)
    bb = last(j)
    idx_i,idx_j = _str_to_idx(param,ii,jj)
    sites_i,sites_j = param.sites[idx_i],param.sites[idx_j]
    _idx_a = findfirst(isequal(aa),sites_i)
    _idx_b = findfirst(isequal(bb),sites_j)
    idx_a = _idx_a === nothing ? 0 : _idx_a
    idx_b = _idx_b === nothing ? 0 : _idx_b
    assoc_view = param[idx_i::Int,idx_j::Int]
    setindex!(assoc_view,val,idx_a::Int,idx_b::Int)
end

function AssocParam(
        name::String,
        components::Vector{String},
        values::MatrixofMatrices,
        allcomponentsites = nothing,
        sourcecsvs = String[],
        sources = String[]
    )
    _values = Compressed4DMatrix(values)
    return AssocParam(name, components, _values, allcomponentsites, sourcecsvs,sources)
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
function Base.convert(::Type{AssocParam{T1}},param::AssocParam{T2}) where {T1<:Number,T2<:Number}
    assoc_values = param.values
    new_assoc_values = T1.(assoc_values.values)
    values = Compressed4DMatrix(new_assoc_values,assoc_values.outer_indices,assoc_values.inner_indices,assoc_values.outer_size,assoc_values.inner_size)
    return AssocParam(param.name,param.components,values,param.sites,param.sourcecsvs,param.sources)
end

function Base.convert(::Type{AssocParam{Bool}},param::AssocParam{<:Union{Int,Float64}})
    assoc_values = param.values
    #@assert all(z->(isone(z) | iszero(z)),assoc_values.values)
    new_assoc_values = Array(Bool.(assoc_values.values))
    values = Compressed4DMatrix(new_assoc_values,assoc_values.outer_indices,assoc_values.inner_indices,assoc_values.outer_size,assoc_values.inner_size)

    return AssocParam(param.name,param.components,values,param.sites,param.sourcecsvs,param.sources)
end

function Base.convert(::Type{AssocParam{String}},param::AssocParam{<:AbstractString})
    assoc_values = param.values
    new_assoc_values = String.(assoc_values.values)
    values = Compressed4DMatrix(new_assoc_values,assoc_values.outer_indices,assoc_values.inner_indices,assoc_values.outer_size,assoc_values.inner_size)
    return AssocParam(param.name,param.components,values,param.sites,param.sourcecsvs,param.sources)
end

#trying to break stack overflow on julia 1.6
function Base.convert(::Type{AssocParam{String}},param::AssocParam{String})
    return param
end