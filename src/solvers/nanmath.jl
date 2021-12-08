@inline log(x) = Base.log(x)

@inline function log(x::T) where T <:Real
        _0 = zero(x)
        ifelse(x>=zero(x),Base.log(max(_0,x)),_0/_0)
    end

@inline sqrt(x) = Base.sqrt(x)
@inline function sqrt(x::T) where T <:Real
    _0 = zero(x)
    ifelse(x>=zero(x),Base.sqrt(max(_0,x)),_0/_0)
end