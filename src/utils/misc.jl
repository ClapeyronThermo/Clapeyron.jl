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
"""
    dnorm(x,y,p)

Equivalent to `norm((xi-yi for (xi, yi) in zip(x, y)), p)`
"""
function dnorm(x,y,p)
    return norm((xi-yi for (xi, yi) in zip(x, y)), p)
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
                "User-Agent" => "https://github.com/ypaul21/Clapeyron.jl"]
    
    url = "https://api.crossref.org/v1/works/" * doi * "/transform"
    out = IOBuffer()
    try
        r = Base.CoreLogging.with_logger(Base.CoreLogging.NullLogger()) do
            Downloads.request(url, output = out, method = "GET",headers = headers)
        end
        if r.status == 200
            res = String(take!(out))
        else  
            res =  ""
        end
        DOI2BIB_CACHE[doi] = res
        return res
    catch err
        @show err
        @warn "no internet connection"
        return ""
    end
end