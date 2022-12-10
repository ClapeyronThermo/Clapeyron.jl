@inline log(x) = Base.log(x)

@inline function log(x::T) where T <:Real
        _0 = zero(x)
        ifelse(x>=_0,Base.log(max(_0,x)),_0/_0)
    end

@inline sqrt(x) = Base.sqrt(x)
@inline function sqrt(x::T) where T <:Real
    _0 = zero(x)
    ifelse(x>=_0,Base.sqrt(max(_0,x)),_0/_0)
end

@inline log1p(x) = Base.log1p(x)

@inline function log1p(x::T) where T <:Real
        _m1 = -one(x)
        _0 = zero(x)
        ifelse(x>=_m1,Base.log1p(max(_m1,x)),_0/_0)
    end

const basepow = Base.:^

@inline function ^(x::Real, y::Real)
    x,y,_nan = promote(x,y,zero(x)/zero(x)) #this will make pow type-stable, at the cost of losing Integer exponentiation
    z = ifelse(x>=zero(x),x,_nan)
    return basepow(z,y)
    end

@inline function ^(x::Real, y::Int)
    return basepow(x,y)
    end
    
^(x,y) = basepow(x,y)
    