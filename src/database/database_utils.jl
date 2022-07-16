const NumberOrString = Union{Union{T1,Missing},Union{T2,Missing}} where {T1 <: AbstractString, T2 <: Number}

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
function getpaths(location::AbstractString; relativetodatabase::Bool=false)
    # We do not use realpath here directly because we want to make the .csv suffix optional.
    filepath = relativetodatabase ? normpath(dirname(pathof(Clapeyron)), "..", "database", location) : location
    isfile(filepath) && return [realpath(filepath)]
    isfile(filepath * ".csv") && return [realpath(filepath * ".csv")]
    if !isdir(filepath)
        relativetodatabase ? error("The path ", location, " does not exist in the Clapeyron database.") :
            error("The path ", location, " does not exist.")
    end
    files = joinpath.(filepath, readdir(filepath))
    return realpath.(files[isfile.(files) .& (getfileextension.(files) .== "csv")])
end

function flattenfilepaths(locations,userlocations)
    res = 
                vcat(
                    reduce(vcat,getpaths.(locations; relativetodatabase=true),init = String[]),
                    reduce(vcat,getpaths.(userlocations),init = String[]),
                    String[]
                    )
    return res
end

function getline(filepath::AbstractString, selectedline::Int)
    # Simple function to return text from filepath at selectedline.
    open(filepath) do file
        linecount = 1
        for line ∈ eachline(file)
            linecount == selectedline && return line
            linecount += 1
        end
        error("Selected line number exceeds number of lines in file")
    end
end

function normalisestring(str, isactivated::Bool=true; tofilter::Regex=r"[ \-\_]", changecase::Bool=true)
    !isactivated && return string(str)
    normalisedstring = replace(str, tofilter => "")
    changecase && (normalisedstring = lowercase(normalisedstring))
    return normalisedstring
end

function _indexin(query,list,separator)
    querydict = Dict(v => k for (k,v) in pairs(query))
    return _indexin(querydict,list,separator,keys(list))
end

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
                if startswith(list_i,ki) #starts with string
                    k_with_separator = ki * separator
                elseif endswith(list_i,ki)
                    k_with_separator = separator * ki
                else
                    k_with_separator = separator * ki * separator  
                end
                if occursin(k_with_separator,list_i)
                    push!(res,k)
                    push!(comp_res,query[ki])   
                end
            end
        end
    end
    return res,comp_res
end

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

function info_color(text)
    colors = Base.text_colors
    red = colors[:bold] * colors[:cyan]
    reset = colors[:normal]
    return red * text * reset
end


