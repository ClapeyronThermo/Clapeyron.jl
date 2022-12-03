

struct ∂Tag end

@inline function derivative(f::F, x::R) where {F,R<:Real}
    return ForwardDiff.derivative(f,x)
end

@inline function gradient2(f::F, x1::R,x2::R) where {F,R<:Real}
    T = typeof(ForwardDiff.Tag(f, R))
    _1 = one(R)
    _0 = zero(R)
    dual1 = ForwardDiff.Dual{T,R,2}(x1, ForwardDiff.Partials((_1,_0)))
    dual2 = ForwardDiff.Dual{T,R,2}(x2, ForwardDiff.Partials((_0,_1)))
    out = f(dual1,dual2)
    ∂out = ForwardDiff.partials(out)
    return SVector(∂out.values)
end

function gradient2(f::F,x1::R1,x2::R2) where{F,R1<:Real,R2<:Real}
    y1,y2 = promote(x1,x2)
    return gradient2(f,y1,y2)
end

"""
    f∂f(f,x)

returns f and ∂f/∂x evaluated in `x`, using `ForwardDiff.jl`, `DiffResults.jl` and `StaticArrays.jl` to calculate everything in one pass.
"""
@inline function f∂f(f::F, x::R) where {F,R<:Real}
    T = typeof(ForwardDiff.Tag(f, R))
    out = f(ForwardDiff.Dual{T,R,1}(x, ForwardDiff.Partials((one(R),))))
    return ForwardDiff.value(out),  ForwardDiff.extract_derivative(T, out)
end

"""
    f∂f∂2f(f,x)

returns f,∂f/∂x,and ∂²f/∂²x and evaluated in `x`, using `ForwardDiff.jl`, `DiffResults.jl` and `StaticArrays.jl` to calculate everything in one pass.
"""
@inline function f∂f∂2f(f::F,x::R) where {F,R<:Real}
    T = typeof(ForwardDiff.Tag(f, R))
    out = ForwardDiff.Dual{T,R,1}(x, ForwardDiff.Partials((one(R),)))
    _f,_df = f∂f(f,out)
    fx = ForwardDiff.value(_f)
    dfx = ForwardDiff.partials(_f).values[1]
    d2fx = ForwardDiff.partials(_df).values[1]
    return (fx,dfx,d2fx)
end

"""
    fgradf2(f,x1,x2)

returns f and ∇f(x),using `ForwardDiff.jl`
"""
function fgradf2(f::F,x1::R1,x2::R2) where{F,R1<:Real,R2<:Real}
    y1,y2 = promote(x1,x2)
    return fgradf2(f,y1,y2)
end

@inline function fgradf2(f::F,x1::R,x2::R) where{F,R<:Real}
    T = typeof(ForwardDiff.Tag(f, R))
    _1 = one(R)
    _0 = zero(R)
    dual1 = ForwardDiff.Dual{T,R,2}(x1, ForwardDiff.Partials((_1,_0)))
    dual2 = ForwardDiff.Dual{T,R,2}(x2, ForwardDiff.Partials((_0,_1)))
    out = f(dual1,dual2)
    ∂out = ForwardDiff.partials(out)
    return ForwardDiff.value(out),SVector(∂out.values)
end

#Manual implementation of an hyperdual.
@inline function ∂2(f::F,x1::R,x2::R) where {F,R<:Real}
    T = typeof(ForwardDiff.Tag(f, R))
    _1 = one(R)
    _0 = zero(R)
    dual1 = ForwardDiff.Dual{T,R,2}(x1, ForwardDiff.Partials((_1,_0)))
    dual2 = ForwardDiff.Dual{T,R,2}(x2, ForwardDiff.Partials((_0,_1)))  
    _f,_df = fgradf2(f,dual1,dual2)
    fx = ForwardDiff.value(_f)
    df1,df2 = _df[1],_df[2]
    df = SVector(df1.value , df2.value)
    d2fdx2, d2fdxdy = df1.partials.values
          _, d2fdy2 = df2.partials.values
    d2f = SMatrix{2}(d2fdx2,d2fdxdy,d2fdxdy,d2fdy2)
    return (fx,df,d2f)
end

function ∂2(f::F,x1::R1,x2::R2) where{F,R1<:Real,R2<:Real}
    y1,y2 = promote(x1,x2)
    return ∂2(f,y1,y2)
end

@inline function ∂J2(f::F,x1::R,x2::R) where {F,R<:Real}
    T = typeof(ForwardDiff.Tag(f, R))
    _1 = one(R)
    _0 = zero(R)
    dual1 = ForwardDiff.Dual{T,R,2}(x1, ForwardDiff.Partials((_1,_0)))
    dual2 = ForwardDiff.Dual{T,R,2}(x2, ForwardDiff.Partials((_0,_1)))  
    _f,_df = fgradf2(f,dual1,dual2)
    fx = ForwardDiff.value(_f)
    df1,df2 = _df[1],_df[2]
    df = SVector(df1.value , df2.value)
    d2fdx2, d2fdxdy = df1.partials.values
          _, d2fdy2 = df2.partials.values
    d2f = SMatrix{2}(d2fdx2,d2fdxdy,d2fdxdy,d2fdy2)
    return (fx,df,d2f)
end

function autochunk(x)
    k = ForwardDiff.pickchunksize(length(x))
    return ForwardDiff.Chunk{k}()
end