const BigNaN = big"NaN"
@inline nan(::T) where T = nan(T)
@inline nan(::Type{Float16}) = NaN16
@inline nan(::Type{Float32}) = NaN32
@inline nan(::Type{Float64}) = NaN64
#TODO: define nan for integer types
@inline nan(::Type{BigFloat}) = BigNaN
@inline nan(::Type{BigInt}) = BigNaN
@inline nan(::Type{T}) where T<:Real = zero(T)/zero(T)
#check if this gives out the correct behaviour.
@inline nan(::Type{T}) where T <: Rational{V} where V = nan(V)

@inline log(x) = Base.log(x)

@inline function log(x::T) where T <:Real
        _0 = zero(T)
        res = Base.log(max(_0,x))
        ifelse(x>=_0,res,nan(res))
    end

@inline sqrt(x) = Base.sqrt(x)
@inline function sqrt(x::T) where T <:Real
    _0 = zero(x)
    res = Base.sqrt(max(_0,x))
    ifelse(x>=_0,res,nan(res))
end

@inline log1p(x) = Base.log1p(x)

@inline function log1p(x::T) where T <:Real
        _m1 = -one(x)
        _0 = zero(x)
        res = Base.log1p(max(_m1,x))
        ifelse(x>=_m1,res,nan(res))
    end

const basepow = Base.:^

@inline function ^(x::Real, y::Real)
    x,y,_nan = promote(x,y,nan(x*y)) #this will make pow type-stable, at the cost of losing Integer exponentiation
    z = ifelse(x>=zero(x),x,_nan)
    return basepow(z,y)
    end

@inline function ^(x::Real, y::Int)
    return basepow(x,y)
    end

^(x,y) = basepow(x,y)
