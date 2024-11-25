#=
This file contains functionalities present in newer versions of julia, but not on LTS,
that are used by this package.
=#
@static if isdefined(Base,Symbol("@assume_effects"))
    macro pure(ex)
        esc(:(Base.@assume_effects :foldable $ex))
    end
else
    macro pure(ex)
        esc(:(Base.@pure $ex))
    end
end

@eval macro $(Symbol("const"))(field)
    if VERSION >= v"1.8.0-DEV.1148"
        Expr(:const, esc(field))
    else
        return esc(field)
    end
end

#Copied from ArrayInterfaceCore.jl
@pure __parameterless_type(T) = Base.typename(T).wrapper

"""
    parameterless_type(x)
given a type or instance of type, it will return the type without any parameters attached.
## Examples:
```julia-repl
julia> Clapeyron.parameterless_type(Vector{Float64})
Array
julia> Clapeyron.parameterless_type(rand(1:5,10))
Array
```
"""
parameterless_type(x) = parameterless_type(typeof(x))
parameterless_type(x::Type) = __parameterless_type(x)

#this is never used in a critical path, so we just use a default copying method
@static if !isdefined(Base,:keepat!)
    function keepat!(a,inds)
        b = a[inds]
        resize!(a,length(b))
        a .= b
        return a
    end
end

@static if !isdefined(Base,:eachsplit)
    eachsplit(str::AbstractString, dlm; limit::Integer=0, keepempty::Bool=true) = split(str,dlm;limit,keepempty)
    eachsplit(str::AbstractString; limit::Integer=0, keepempty::Bool=false)  = split(str;limit,keepempty)
end

split_2(str) = NTuple{2}(eachsplit(str, limit=2))
split_2(str,dlm) = NTuple{2}(eachsplit(str,dlm, limit=2))




function _vecparser_eltype(vals)
    for val in eachsplit(vals,' ')
        if isnothing(tryparse(Int,val))
            return Float64 
        end
    end
    return Int
end

function _vecparser(T::Type{X},vals::String,dlm = ' ') where X <: Union{Int,Float64}
    strip_vals = strip(vals,('[',']'))
    res = Vector{T}(undef,0)
    for strval in eachsplit(strip_vals,dlm,keepempty = false)
        val = tryparse(T,strval)
        if !isnothing(val)
            push!(res,val)
        else
            colors = Base.text_colors
            red = colors[:bold] * colors[:red]
            reset = colors[:normal]
            errval = red * strval * reset
            error("cannot parse $errval as a number in $vals")
        end
    end
    return res
end

function _vecparser(vals::String,dlm = ' ')
    strip_vals = strip(vals,('[',']'))
    T = _vecparser_eltype(strip_vals)
    return _vecparser(T,vals,dlm)
end
