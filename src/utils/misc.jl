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

function doi2bib(doi::String)
    headers = ["Accept" => "application/x-bibtex"]
    url = "https://dx.doi.org/" * doi
    @show url
    out = IOBuffer()
    r = Downloads.request(url, output = out, method = "GET",headers = headers)
    if r.status == 200
        res = String(take!(out))
    else
        return ""
    end
end