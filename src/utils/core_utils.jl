#Copied from ArrayInterfaceCore.jl
Base.@assume_effects :foldable __parameterless_type(T) = Base.typename(T).wrapper

"""
    parameterless_type(x)
Given a type or instance of type, it will return the type without any parameters attached.
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

split_2(str) = NTuple{2}(eachsplit(str, limit=2))
split_2(str,dlm) = NTuple{2}(eachsplit(str,dlm, limit=2))

function fma_evalpoly(x::T1,pol::NTuple{N,T2}) where {T1,N,T2}
    fx = fma(pol[end],x,pol[end - 1])
    @inbounds for i in 2:(N-1)
        fx = fma(fx,x,pol[end-i])
    end
    return fx
end

function show_pairs(io,keys,vals=nothing,separator="",f_print = print;quote_string = true,pair_separator = '\n',prekey = ifelse(pair_separator === '\n'," ",""))
    if length(keys) == 0
        return nothing
    end
    if vals === nothing #useful for printing only keys
        vals = Iterators.repeated("")
    end
    i = 0
    for (k,v) in zip(keys,vals)
        i += 1
        if i > 1
            print(io,pair_separator)
        end
        if quote_string
            quot = '\"'
            print(io,prekey,quot,k,quot,separator)
        else
            print(io,prekey,k,separator)
        end
        f_print(io,v)
    end
end

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

show_default(io::IO,arg) = Base.show_default(io,arg)
show_default(io::IO,mime::MIME"text/plain",arg) = Base.show_default(io,arg)

function show_as_namedtuple(io::IO,x)
    compact_io = IOContext(io, :compact => true)
    print(io,parameterless_type(x))
    print(io,"(")
    names = fieldnames(typeof(x))
    l = length(names)
    equal = " = "
    comma = ", "
    for i in 1:l
        print(io,names[i])
        print(io,equal)
        print(compact_io,getfield(x,i))
        if i != l
            print(io,comma) 
        end
    end
    print(io,")")
end

_zero(t::Number) = zero(t)
_zero(x::T)  where T = _zero(T)
_zero(::Type{T}) where T <: Number = zero(T)
_zero(::Type{String}) = ""
_zero(::Type{Missing}) = missing
_zero(::Type{T}) where T <: AbstractString = T("")
_zero(::Type{T}) where T <:Union{T1,Missing} where T1 = _zero(nonmissingtype(T))


_iszero(t::Number) = iszero(t)
_iszero(::Missing) = true
_iszero(t::AbstractString) = isempty(t)

function raw_values end
raw_values(x) = x

"""
    param_from_values(newval,param)
    param_from_values(f,param)

Given a Clapeyron parameter `param`, returns another Clapeyron parameter with their raw values replaced with `newval`.
The modification can occur inplace, or allocate a different value structure.
Also works in function form:

```julia
x = SingleParam("bb",["a","b"]) #singleparam filled with zeros
Clapeyron.param_from_values(x) do values
    values .+= 1
end
```

For nested parameters it work until the flat storage. (the storage in `Compressed4DMatrix`, for example)

"""
param_from_values(f::C,param::P) where {C <: Base.Callable,P} = _param_from_values(f(raw_values(param)),param)
param_from_values(x,param::P) where P = _param_from_values(x,param)
#=
"""
    concrete(x)

Given an array of heterogeneous values, return an array of concrete values.
"""
concrete(x::Vector{Float64}) = x
concrete(x::Vector{Int64}) = x
concrete(x::Vector{String}) = x
concrete(x::Vector{Bool}) = x
concrete(x) = convert(AbstractArray{mapreduce(typeof, promote_type, x)}, x)=#
