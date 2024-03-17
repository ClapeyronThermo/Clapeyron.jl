"""
    vec2(x1,x2,opt=true)
Generates a correct 2-length static array `[x1,x2]`, with support for non-isbits types
"""
function vec2(x1,x2,opt = true)
    V01,V02,_ = promote(x1,x2,opt)
    if V01 isa Base.IEEEFloat # MVector does not work on non bits types, like BigFloat
        return MVector((V01,V02))
    else
        return SizedVector{2,typeof(V01)}((V01,V02))
    end
end

function vec3(x1,x2,x3,opt = true)
    V01,V02,V03,_ = promote(x1,x2,x3,opt)
    if V01 isa Base.IEEEFloat # MVector does not work on non bits types, like BigFloat
        return MVector(V01,V02,V03)
    else
        return SizedVector{2,typeof(V01)}((V01,V02,V03))
    end
end

function svec2(x1,x2,opt = true)
    V01,V02,_ = promote(x1,x2,opt)
    return SVector{2}(V01,V02)
end

function svec3(x1,x2,x3,opt = true)
    V01,V02,V03,_ = promote(x1,x2,x3,opt)
    return SVector{3}(V01,V02,V03)
end

"""
    _parse_kv(str,dlm)

Parses a key-value pair from a string, returns the key and the value:

## Example
```
julia> Clapeyron._parse_kv("  a = b  ",'=')
("a", "b")
```

"""
function _parse_kv(str,dlm)
    _k,_v = split_2(str,dlm)
    return strip(_k),strip(_v)
end

DOI2BIB_CACHE = Dict{String,String}()

"""
    doi2bib(doi::String)

given a DOI identifier, returns a BibTeX entry. requires an internet connection. if the value is not found, returns an empty string. results are cached in `Clapeyron.DOI2BIB_CACHE`

## Example

```julia-repl
julia> Clapeyron.doi2bib("10.1063/1.5136079")
"@article{Aasen_2020,\n\tdoi = {10.1063/1.5136079},\n\turl = {https://doi.org/10.1063%2F1.5136079},\n\tyear = 2020,\n\tmonth = {feb},\n\tpublisher = {{AIP} Publishing},\n\tvolume = {152},\n\tnumber = {7},\n\tpages = {074507},\n\tauthor = {Ailo Aasen and Morten Hammer and Erich A. Müller and {\\O}ivind Wilhelmsen},\n\ttitle = {Equation of state and force fields for Feynman{\\textendash}Hibbs-corrected Mie fluids. {II}. Application to mixtures of helium, neon, hydrogen, and deuterium},\n\tjournal = {The Journal of Chemical Physics}\n}"
```


"""
function doi2bib(doi::String)

    if haskey(DOI2BIB_CACHE,doi)
        return DOI2BIB_CACHE[doi] #caching requests
    end

    headers = ["Accept"=>"application/x-bibtex",
                "charset" => "utf-8",
                "User-Agent" => "https://github.com/ClapeyronThermo/Clapeyron.jl"]

    url = "https://api.crossref.org/v1/works/" * doi * "/transform"
    out = IOBuffer()
    try
        r = Base.CoreLogging.with_logger(Base.CoreLogging.NullLogger()) do
            Downloads.request(url, output = out, method = "GET",headers = headers)
        end
        if r.status == 200
            res = strip(String(take!(out)))
        else
            res = ""
        end
        DOI2BIB_CACHE[doi] = res
        return res
    catch err
        @show err
        @warn "no internet connection"
        return ""
    end
end

"""
    evalexppoly(x,n,v)

Returns ∑nᵢx^vᵢ.
"""
function evalexppoly(x,n,v)
    res = zero(x*first(n)*first(v))
    logx = log(x)
    for i in 1:length(n)
        res += n[i]*exp(logx*v[i])
    end
    return res
end

function cached_indexin(a, b, bdict)
    inds = keys(b)
    #bdict = Dict{eltype(b),eltype(inds)}()
    empty!(bdict)
    for (val, ind) in zip(b, inds)
        get!(bdict, val, ind)
    end
    return Union{eltype(inds), Nothing}[
        get(bdict, i, nothing) for i in a
    ]
end

format_components(str::String) = [str]
format_components(str::Tuple) = format_components(first(str))
format_components(str::Pair) = format_components(first(str))
format_components(str::AbstractString) = format_components(String(str))
format_components(str::Vector{String}) = str
format_components(str) = map(format_component_i,str)
format_component_i(str::AbstractString) = String(str)
format_component_i(str::String) = str
format_component_i(x::Tuple) = first(x)
format_component_i(x::Pair) = first(x)

format_gccomponents(str::Tuple) = [str]
format_gccomponents(str::Pair) = [str]
format_gccomponents(str) = str
format_gccomponents(str::String) = [str]
format_gccomponents(str::AbstractString) = format_components(String(str))
format_gccomponents(str::Vector{String}) = str
#used by MultiComponentFlash extension