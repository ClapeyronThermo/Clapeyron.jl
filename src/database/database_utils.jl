const NumberOrString = Union{Union{T1,Missing},Union{T2,Missing}} where {T1 <: AbstractString, T2 <: Number}

const DB_PATH = normpath(Base.pkgdir(Clapeyron),"database")

const SHORT_PATHS = Dict{String,String}(
    "DB" => DB_PATH
)

const SPECIAL_IDENTIFIERS = ["@REPLACE"]

const SKIP_GETPATHS =  ("Clapeyron Database File", #a raw CSV file
                        "Clapeyron Estimator")

"""
    getfileextension(filepath)
A quick helper to get the file extension of any given path (without the dot).
# Examples
```julia-repl
julia> getfileextension("~/Desktop/text.txt")
"txt"
```
"""
function getfileextension(filepath::AbstractString)
    path,ext = splitext(filepath)
    return String(chop(ext,head=1,tail=0))
end


"""
    getpaths(location; relativetodatabase=false)
Returns database paths that is optionally relative to Clapeyron.jl directory.
If path is a file, then return an Array containing a single path to that file.
If path is a directory, then return an Array containing paths to all csv files in that directory.
# Examples
```julia-repl
julia> getpaths("SAFT/PCSAFT"; relativetodatabase=true)
3-element Array{String,1}:
 "/home/user/.julia/packages/Clapeyron.jl/xxxxx/database/SAFT/PCSAFT/data_PCSAFT_assoc.csv"
 "/home/user/.julia/packages/Clapeyron.jl/xxxxx/database/SAFT/PCSAFT/data_PCSAFT_like.csv"
 "/home/user/.julia/packages/Clapeyron.jl/xxxxx/database/SAFT/PCSAFT/data_PCSAFT_unlike.csv"
```
"""
function getpaths(location::AbstractString; relativetodatabase::Bool=false)::Vector{String}
    # We do not use realpath here directly because we want to make the .csv suffix optional.
    is_inline_csv(location) && return [location]
    if startswith(location,"@REPLACE")
        filepath = chop(location,head = 9, tail = 0)
        result = getpaths(filepath)
        rr = ["@REPLACE" * Base.Filesystem.path_separator * res for res in result]
        return rr
    end
    if relativetodatabase
        new_loc = normpath("@DB",location) #we suppose that the database is never at the root of a windows drive
    else
        new_loc = location
    end
    return _getpaths(new_loc)
end

function _getpaths(location,special_parse = true)
    location == "@REMOVEDEFAULTS" && return [location]
    if special_parse && startswith(location,'@')
        locs = splitpath(location)
        first_identifier = locs[1]
        if startswith(first_identifier,'@')
            raw_first_identifier = chop(first_identifier,head = 1,tail = 0)
            if haskey(SHORT_PATHS,raw_first_identifier)
                locs[1] = SHORT_PATHS[raw_first_identifier]
                return _getpaths(join(locs,Base.Filesystem.path_separator))
            else
                return _getpaths(location,false)
            end
        end
    end
    filepath = location
    isfile(filepath) && return [realpath(filepath)]
    
    #if we want to parse jsons, this is ambiguous.
    #isfile(filepath * ".csv") && return [realpath(filepath * ".csv")]
    
    #=
    this should fail at the CSV reader stage
    if !isdir(filepath)
        relativetodatabase ? error("The path ", location, " does not exist in the Clapeyron database.") :
            error("The path ", location, " does not exist.")
    end =#
    files = readdir(filepath,join = true) #this returns the full (non-normalized) path
    filter!(isfile,files) #remove folders, the reader is not recursive
    filter!(f -> getfileextension(f) in ("csv","json"),files)
    map!(realpath,files,files)
    return files
end

flattenfilepaths(locations) = flattenfilepaths(locations,String[])

function flattenfilepaths(locations,userlocations::Vector{String})
    if length(locations) == 0 && length(userlocations) == 0
        return String[]
    end
    defaultpaths = reduce(vcat,getpaths.(locations; relativetodatabase=true),init = String[])
    userpaths = reduce(vcat,getpaths.(userlocations),init = String[])
    idx = findfirst(isequal("@REMOVEDEFAULTS"),userpaths)
    if !isnothing(idx)
        defaultpaths = String[]
        popat!(userpaths,idx)
    end
    return vcat(defaultpaths,userpaths,String[])
end

flattenfilepaths(locations,userlocations::AbstractString) = flattenfilepaths(locations,[userlocations])

getpath(location;relativetodatabase = true) = only(getpaths(location; relativetodatabase))

Base.@nospecialize
function flattenfilepaths(locations,userlocations)
    return String[]
end
Base.@specialize

function getline(filepath::AbstractString, selectedline::Int)
    is_inline_csv(filepath) && return getline(IOBuffer(filepath),selectedline)
    open(filepath) do file
       _getline(file,selectedline)
    end
end

getline(file::IO,selectedline::Int) = _getline(file,selectedline)

function _getline(file, selectedline::Int)
    # Simple function to return text from filepath at selectedline.
    linecount = 1
    for line ∈ eachline(file)
        linecount == selectedline && return line
        linecount += 1
    end
    error("Selected line number exceeds number of lines in file")
end

function normalisestring(str, isactivated::Bool=true; tofilter = ' ')
    ismissing(str) && return ""
    if !isactivated
        str isa String && return str::String
        return string(str)::String
    end
    res = Base.Unicode.normalize(str,casefold=true,stripmark=true)
    return replace(res, tofilter => "")
end

function is_inline_csv(filepath)
    return any(startswith(filepath,kw) for kw in SKIP_GETPATHS)
end

function _indexin(query,list,separator)
    querydict = Dict(v => k for (k,v) in pairs(query))
    return _indexin(querydict,list,separator,keys(list))
end

if isdefined(Base,:eachsplit)
function _indexin(query,list,separator,indices)
    kq = keys(query)
    res = zeros(Int,0)
    comp_res = zeros(Int,0)
    sizehint!(res,2*length(kq))
    sizehint!(comp_res,2*length(kq))
    for k in indices
        list_i = list[k]
        for val in eachsplit(list_i,separator,keepempty = false)
            if val in kq
                push!(res,k)
                push!(comp_res,query[val])
            end
        end
    end
    return res,comp_res
end

else
    function _indexin(query,list,separator,indices)
        kq = keys(query)
        res = zeros(Int,0)
        comp_res = zeros(Int,0)
        sizehint!(res,2*length(kq))
        sizehint!(comp_res,2*length(kq))
        for k in indices
            list_i = list[k]
            if !occursin(separator,list_i) #simple format
                if list_i in kq
                    push!(res,k)
                    push!(comp_res,query[list_i])
                end
            else #separator format
                for ki in kq
                    if startswith(list_i,ki * separator) #starts with string
                        push!(res,k)
                        push!(comp_res,query[ki])
                    elseif endswith(list_i,separator * ki)  #ends with string
                        push!(res,k)
                        push!(comp_res,query[ki])
                    elseif occursin(separator * ki * separator,list_i) #string in between
                        push!(res,k)
                        push!(comp_res,query[ki])
                    end
                end
            end
        end
        return res,comp_res
    end
end
#=
function _indexin(query,list,separator,indices)
    kq = keys(query)
    res = zeros(Int,0)
    comp_res = zeros(Int,0)
    sizehint!(res,2*length(kq))
    sizehint!(comp_res,2*length(kq))
    for k in indices
        list_i = list[k]
        if !occursin(separator,list_i) #simple format
            if list_i in kq
                push!(res,k)
                push!(comp_res,query[list_i])
            end
        else #separator format
            for ki in kq
                if startswith(list_i,ki * separator) #starts with string
                    push!(res,k)
                    push!(comp_res,query[ki])
                elseif endswith(list_i,separator * ki)  #ends with string
                    push!(res,k)
                    push!(comp_res,query[ki])
                elseif occursin(separator * ki * separator,list_i) #string in between
                    push!(res,k)
                    push!(comp_res,query[ki])
                end
            end
        end
    end
    return res,comp_res
end =#
function defaultmissing(array::Array{<:Number},defaultvalue = zero(eltype(array)))
    return deepcopy(array),Array(ismissing.(array))
end

function defaultmissing(array::Array{<:AbstractString},defaultvalue = "")
    return string.(array),Array(ismissing.(array))
end

function defaultmissing(array::Array{String},defaultvalue = "")
    return deepcopy(array),Array(ismissing.(array))
end

function defaultmissing(array::Array{Union{Missing, T}},defaultvalue="") where T <:AbstractString
    return string.(coalesce.(array,defaultvalue)),Array(ismissing.(array))
end

function defaultmissing(array::Array{Union{Missing, T}},defaultvalue=zero(eltype(array))) where T<:Number
    return coalesce.(array,defaultvalue),Array(ismissing.(array))
end

function defaultmissing(array::Array{Union{Missing, Bool}},defaultvalue=zero(eltype(array)))
    return [coalesce(i,defaultvalue) for i in array],Array(ismissing.(array))
end

#if an array with only missings is passed, the Resulting ClapeyronParam will be
#of the type that this function returns
function defaultmissing(array::Array{Missing},defaultvalue=0.0)
    return coalesce.(array,defaultvalue),Array(ismissing.(array))
end

function defaultmissing(array::Array{Union{Missing,T1,T2}},defaultvalue="") where {T1 <:Number,T2<:AbstractString}
    return string.(coalesce.(array,defaultvalue)),Array(ismissing.(array))
end

function defaultmissing(array::Array{Any},defaultvalue="")
    return string.(coalesce.(array,defaultvalue)),Array(ismissing.(array))
end

function defaultmissing(array::Array{T},defaultvalue::T2) where T<:Union{T2,Missing} where T2
    coalesce.(array,Ref(defaultvalue)),Array(ismissing.(array))
end
function defaultmissing(array)
    throw("Unsupported array element type  $(typeof(array))")
end

_zero(t::Number) = zero(t)
_zero(::String) = ""
_zero(::Missing) = missing
function _zero(x::Type{T})  where T<:Number
    return zero(T)
end
_zero(::Type{String}) = ""
_zero(::Type{Missing}) = missing
function _zero(::Type{T}) where T <:Union{T1,Missing} where T1
    return missing
end

_iszero(t::Number) = iszero(t)
_iszero(::Missing) = true
_iszero(t::AbstractString) = isempty(t)

"""
    singletopair(params::Vector,outputmissing=zero(T))
Generates a square matrix, filled with "zeros" (considering the "zero" of a string, a empty string).
The generated matrix will have the values of `params` in the diagonal.
If missing is passed, the matrix will be filled with `missing`
"""
function singletopair(params::Vector{T1},::T2 =_zero(T1)) where {T1,T2}
    len = length(params)
    T = Union{T1,T2}
    output = Matrix{T}(undef,len,len)
    fill!(output,_zero(T))
    @inbounds  for i in 1:len
        output[i,i] = params[i]
    end
    return output
end

##Error and info display utils
function error_color(text)
    colors = Base.text_colors
    red = colors[:bold] * colors[:red]
    reset = colors[:normal]
    return red * text * reset
end

error_color(symbol::Symbol) = error_color(":" * string(symbol))

function info_color(text)
    colors = Base.text_colors
    red = colors[:bold] * colors[:cyan]
    reset = colors[:normal]
    return red * text * reset
end

info_color(symbol::Symbol) = info_color(":" * string(symbol))


function userlocation_merge(loc1,loc2)
    if isempty(loc2)
        return loc1
    elseif isempty(loc1)
        return loc2
    elseif loc1 isa Vector{String} && loc2 isa Vector{String}
        return append!(loc1,loc2)
    elseif loc1 isa Vector{String} && length(loc1) == 0
        return loc2
    elseif loc1 isa Vector{String} && can_nt(loc2)
        return loc2
    elseif can_nt(loc1) && can_nt(loc2)
        loc1p = to_nt(loc1)
        loc2p = to_nt(loc2)
        loc0 = Dict(pairs(loc1p))
        for k in keys(loc2p)
            loc0[k] = loc2p[k]
        end
        return loc0
    else
        throw(ArgumentError("invalid userlocations combination: old: $loc1, new: $loc2"))
    end
end

critical_data() = ["properties/critical.csv"]
mw_data() = ["properties/molarmass.csv"]

function by_cas(caslist)
    cas = format_components(caslist)
    params = getparams(cas,["properties/identifiers.csv"],species_columnreference = "CAS",ignore_headers = String[],ignore_missing_singleparams = String["SMILES","inchikey","species"])
    species = params["species"].values
    @show
    for (i,sp) in pairs(species)
        if occursin("~|~",sp)
            x,_ = eachsplit(sp,"~|~")
            species[i] = x
        end
    end
    return species
end

function cas(components)
    components = format_components(components)
    params = getparams(components,["properties/identifiers.csv"],ignore_headers = String["SMILES"],ignore_missing_singleparams = ["CAS"])
    return params["CAS"].values
end

function SMILES(components)
    components = format_components(components)
    params = getparams(components,["properties/identifiers.csv"],ignore_headers = String["CAS"])
    return params["SMILES"].values
end