

struct ∂Tag end

@inline function derivative(f::F, x::R) where {F,R<:Real}
    return ForwardDiff.derivative(f,x)
end

@inline function gradient(f::F, x) where {F}
    return ForwardDiff.gradient(f,x)
end

@inline function gradient!(fx::R,f::F, x) where {R,F}
    return ForwardDiff.gradient!(fx,f,x)::R
end

@inline function hessian(f::F, x) where {F}
    return ForwardDiff.hessian(f,x)
end

@inline function gradient2(f::F, x1::R,x2::R) where {F,R<:Real}
    T = typeof(ForwardDiff.Tag(f, R))
    _1 = oneunit(R)
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
    out = f(ForwardDiff.Dual{T,R,1}(x, ForwardDiff.Partials((oneunit(R),))))
    return ForwardDiff.value(out),  ForwardDiff.extract_derivative(T, out)
end

f∂f(f::F) where F = Base.Fix1(f∂f,f)

"""
    f∂f∂2f(f,x)

returns f,∂f/∂x,and ∂²f/∂²x and evaluated in `x`, using `ForwardDiff.jl`, `DiffResults.jl` and `StaticArrays.jl` to calculate everything in one pass.
"""
@inline function f∂f∂2f(f::F,x::R) where {F,R<:Real}
    T = typeof(ForwardDiff.Tag(f, R))
    out = ForwardDiff.Dual{T,R,1}(x, ForwardDiff.Partials((oneunit(R),)))
    _f,_df = f∂f(f,out)
    fx = ForwardDiff.value(_f)
    dfx = ForwardDiff.partials(_f).values[1]
    d2fx = ForwardDiff.partials(_df).values[1]
    return (fx,dfx,d2fx)
end

f∂f∂2f(f::F) where F = Base.Fix1(f∂f∂2f,f)

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
    _1 = oneunit(R)
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
    _1 = oneunit(R)
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

#Manual implementation of an hyperdual.
@inline function J2(f::F,x::SVector{2,R}) where {F,R<:Real}
    T = typeof(ForwardDiff.Tag(f, R))
    _1 = oneunit(R)
    _0 = zero(R)
    x1,x2 = x
    dual1 = ForwardDiff.Dual{T,R,2}(x1, ForwardDiff.Partials((_1,_0)))
    dual2 = ForwardDiff.Dual{T,R,2}(x2, ForwardDiff.Partials((_0,_1)))  
    dx = SVector(dual1,dual2)
    f̄ = f(dx)
    f̄1,f̄2 = f̄[1],f̄[2]
    F̄ = SVector(f̄1.value , f̄2.value)
    df1dx1, df1dx2 = f̄1.partials.values
    df2dx1, df2dx2 = f̄2.partials.values
    J = SMatrix{2}(df1dx1,df2dx1,df1dx2,df2dx2)
    return F̄,J
end

@inline function J3(f::F,x::SVector{3,R}) where {F,R<:Real}
    T = typeof(ForwardDiff.Tag(f, R))
    _1 = oneunit(R)
    _0 = zero(R)
    x1,x2 = x
    dual1 = ForwardDiff.Dual{T,R,3}(x1, ForwardDiff.Partials((_1,_0,_0)))
    dual2 = ForwardDiff.Dual{T,R,3}(x2, ForwardDiff.Partials((_0,_1,_0)))
    dual3 = ForwardDiff.Dual{T,R,3}(x2, ForwardDiff.Partials((_0,_0,_1)))  
    dx = SVector(dual1,dual2,dual3)
    f̄ = f(dx)
    f̄1,f̄2,f̄3 = f̄[1],f̄[2],f̄[3]
    F̄ = SVector(f̄1.value , f̄2.value, f̄3.value)
    df1dx1, df1dx2, df1dx3 = f̄1.partials.values
    df2dx1, df2dx2, df2dx3 = f̄2.partials.values
    df3dx1, df3dx2, df3dx3 = f̄3.partials.values
    J = SMatrix{3}(df1dx1,df2dx1,df3dx1,df1dx2,df2dx2,df3dx2,df1dx3,df2dx3,df3dx3)
    return F̄,J
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

"""
    primalval(x::Real)

returns the primal value of a value. strips all duals from `ForwardDiff`. useful in debugging:

## Example
```julia
function f(x)
    x1 = primalval(x[1])
    @show x1
    @show x[1]
    sum(x)*exp(x[1]) + log(x[end])
end

julia> ForwardDiff.hessian(f,[1.0,2.0,3.0])
x1 = 1.0
x[1] = Dual{ForwardDiff.Tag{typeof(f), Float64}}(Dual{ForwardDiff.Tag{typeof(f), Float64}}(1.0,1.0,0.0,0.0),Dual{ForwardDiff.Tag{typeof(f), Float64}}(1.0,0.0,0.0,0.0),Dual{ForwardDiff.Tag{typeof(f), Float64}}(0.0,0.0,0.0,0.0),Dual{ForwardDiff.Tag{typeof(f), Float64}}(0.0,0.0,0.0,0.0))
3×3 Matrix{Float64}:
 21.7463   2.71828   2.71828
  2.71828  0.0       0.0
  2.71828  0.0      -0.111111
```

"""
#fallback
primalval(x) = x

#scalar
primalval(x::ForwardDiff.Dual) = primalval(ForwardDiff.value(x))

#primaltype(::Type{T}) where T = T
#primaltype(::Type{<:ForwardDiff.Dual{T,R}}) where {T,R} = primaltype(R)

#arrays overload
function primalval(x::AbstractArray{T}) where T <: ForwardDiff.Dual
    return primalval.(x)
end

#=
gradient at index i

=#

struct GradᵢVector{T,V} <: AbstractVector{T}
    i::Int
    val::T
    vector::V
end

function Base.getindex(x::GradᵢVector{T,V},i) where {T,V}
    idx = x.i
    if idx == i
        return x.val
    else
        return convert(T,x.vector[i])
    end
end

Base.length(x::GradᵢVector) = length(x.vector)
Base.size(x::GradᵢVector) = size(x.vector)

function grad_at_i(f::F,x::X,i,TT = eltype(x)) where {F,X <: AbstractVector{R}} where R
    T = typeof(ForwardDiff.Tag(f, TT))
    xᵢ = TT(x[i])
    ∂xᵢ = ForwardDiff.Dual{T}(xᵢ, oneunit(xᵢ))
    ∂x = GradᵢVector(i,∂xᵢ,x)
    fx = f(∂x)
    return ForwardDiff.extract_derivative(T, fx)
end


