const ForwardDiffStatic = Base.get_extension(ForwardDiff,:ForwardDiffStaticArraysExt)
recursive_fd_value(x::Number) = ForwardDiff.value(x)
recursive_fd_value(x::Tuple) = recursive_fd_value.(x)
recursive_fd_value(x::AbstractArray) = recursive_fd_value.(x)

recursive_fd_extract_derivative(T::TT,x::Number) where TT = ForwardDiff.extract_derivative(T,x)
recursive_fd_extract_derivative(T::TT,x::Tuple) where TT = recursive_fd_extract_derivative.(T,x)
recursive_fd_extract_derivative(T::TT,x::Tuple{}) where TT = x

recursive_fd_extract_derivative(T::TT,x::AbstractArray) where TT = recursive_fd_extract_derivative.(T,x)

@inline function derivative(f::F, x::R, TAG = f) where {F,R<:Real}
    T = typeof(ForwardDiff.Tag(TAG, R))
    return recursive_fd_extract_derivative(T, f(ForwardDiff.Dual{T}(x, one(x))))
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

@inline function gradient2(f::F, x1::R,x2::R,tag = f) where {F,R<:Real}
    T = typeof(ForwardDiff.Tag(tag, R))
    _1 = oneunit(R)
    _0 = zero(R)
    dual1 = ForwardDiff.Dual{T,R,2}(x1, ForwardDiff.Partials((_1,_0)))
    dual2 = ForwardDiff.Dual{T,R,2}(x2, ForwardDiff.Partials((_0,_1)))
    out = f(dual1,dual2)
    ∂out = ForwardDiff.partials(out)
    return SVector(∂out.values)
end

function gradient2(f::F,x1::R1,x2::R2,tag = f) where{F,R1<:Real,R2<:Real}
    y1,y2 = promote(x1,x2)
    return gradient2(f,y1,y2,tag)
end

"""
    f∂f(f,x)

Returns f and ∂f/∂x evaluated in `x`, using `ForwardDiff.jl`, `DiffResults.jl` and `StaticArrays.jl` to calculate everything in one pass.
"""
@inline function f∂f(f::F, x::R,tag = f) where {F,R<:Real}
    T = typeof(ForwardDiff.Tag(tag, R))
    out = f(ForwardDiff.Dual{T,R,1}(x, ForwardDiff.Partials((oneunit(R),))))
    return recursive_fd_value(out),  recursive_fd_extract_derivative(T, out)
end


f∂f(f::F) where F = Base.Fix1(f∂f,f)

"""
    f∂f∂2f(f,x)

Returns f,∂f/∂x,and ∂²f/∂²x and evaluated in `x`, using `ForwardDiff.jl`, `DiffResults.jl` and `StaticArrays.jl` to calculate everything in one pass.
"""
@inline function f∂f∂2f(f::F,x::R,tag = f) where {F,R<:Real}
    T = typeof(ForwardDiff.Tag(tag, R))
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

Returns f and ∇f(x), using `ForwardDiff.jl`.
"""
function fgradf2(f::F,x1::R1,x2::R2,tag = f) where{F,R1<:Real,R2<:Real}
    y1,y2 = promote(x1,x2)
    return fgradf2(f,y1,y2,tag)
end

@inline function fgradf2(f::F,x1::R,x2::R,tag = f) where{F,R<:Real}
    T = typeof(ForwardDiff.Tag(tag, R))
    _1 = oneunit(R)
    _0 = zero(R)
    dual1 = ForwardDiff.Dual{T,R,2}(x1, ForwardDiff.Partials((_1,_0)))
    dual2 = ForwardDiff.Dual{T,R,2}(x2, ForwardDiff.Partials((_0,_1)))
    out = f(dual1,dual2)
    ∂out = ForwardDiff.partials(out)
    return ForwardDiff.value(out),SVector(∂out.values)
end

function ∂2(f::F,x1::R1,x2::R2,tag = f) where{F,R1<:Real,R2<:Real}
    y1,y2 = promote(x1,x2)
    return ∂2(f,y1,y2,tag)
end

#Manual implementation of an hyperdual.
@inline function ∂2(f::F,x1::R,x2::R,tag = f) where {F,R<:Real}
    T = typeof(ForwardDiff.Tag(tag, R))
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
@inline function J2(f::F,x::SVector{2,R},tag = f) where {F,R<:Real}
    T = typeof(ForwardDiff.Tag(tag, R))
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

@inline function J3(f::FF,x::SVector{3,R},tag = f) where {FF,R<:Real}
    T = typeof(ForwardDiff.Tag(tag, R))
    _1 = oneunit(R)
    _0 = zero(R)
    x1,x2,x3 = x
    dx1 = ForwardDiff.Dual{T,R,3}(x1, ForwardDiff.Partials((_1,_0,_0)))
    dx2 = ForwardDiff.Dual{T,R,3}(x2, ForwardDiff.Partials((_0,_1,_0)))
    dx3 = ForwardDiff.Dual{T,R,3}(x3, ForwardDiff.Partials((_0,_0,_1)))  
    dx = SVector(dx1,dx2,dx3)
    f̄ = f(dx)
    f̄1,f̄2,f̄3 = f̄[1],f̄[2],f̄[3]
    Fx = SVector(f̄1.value, f̄2.value, f̄3.value)
    df1dx1, df1dx2, df1dx3 = f̄1.partials.values
    df2dx1, df2dx2, df2dx3 = f̄2.partials.values
    df3dx1, df3dx2, df3dx3 = f̄3.partials.values
    Jx = SMatrix{3}(df1dx1,df2dx1,df3dx1,df1dx2,df2dx2,df3dx2,df1dx3,df2dx3,df3dx3)
    return Fx,Jx
end

function FJ_ad(f::F,x::SVector{3,R},tag = f) where {F,R<:Real}
    return J3(f,x,tag)
end

function FJ_ad(f::F,x::SVector{2,R},tag = f) where {F,R<:Real}
    return J2(f,x,tag)
end

function FJ_ad(f::F,x::X,tag = f) where {F,X}
    Jresult = DiffResults.JacobianResult(x)
    result = ForwardDiff.jacobian!(Jresult,f,x)
    Fx = DiffResults.value(result)
    Jx = DiffResults.jacobian(result)
    return Fx,Jx
end


@inline function ∂J2(f::F,x1::R,x2::R,tag = f) where {F,R<:Real}
    T = typeof(ForwardDiff.Tag(tag, R))
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


#trying to fix ForwardDiff$720

@generated function _extract_jacobian(::Type{T}, ydual::ForwardDiff.Partials{M}, x::S) where {T,M,S<:StaticArray}
    N = length(x)
    result = Expr(:tuple, [:(ForwardDiff.partials(T, ydual[$i], $j)) for i in 1:M, j in 1:N]...)
    return quote
        $(Expr(:meta, :inline))
        V = StaticArrays.similar_type(S, ForwardDiff.valtype(eltype($ydual)), Size($M, $N))
        return V($result)
    end
end

@generated function _extract_gradient(::Type{T}, y::Real, x::S) where {T,S<:StaticArray}
    result = Expr(:tuple, [:(ForwardDiff.partials(T, y, $i)) for i in 1:length(x)]...)
    return quote
        $(Expr(:meta, :inline))
        V = StaticArrays.similar_type(S, ForwardDiff.valtype($y))
        return V($result)
    end
end

function static_fgh(result::DiffResults.ImmutableDiffResult, f::F, x::SVector) where {F}
    T = typeof(ForwardDiff.Tag(f, eltype(x)))
    d1 = ForwardDiffStatic.dualize(T, x)
    d2 = ForwardDiffStatic.dualize(T, d1)
    fd2 = f(d2)
    val = ForwardDiff.value(T,ForwardDiff.value(T,fd2))
    grad = _extract_gradient(T,ForwardDiff.value(T,fd2), x)
    hess = _extract_jacobian(T,ForwardDiff.partials(T,fd2), x)
    result = DiffResults.hessian!(result, hess)
    result = DiffResults.gradient!(result, grad)
    result = DiffResults.value!(result, val)
    return result
end

#obtaining GradientConfig from HessianConfig

function _GradientConfig(hconfig::ForwardDiff.HessianConfig{T,V,N}) where {T,V,N}
    gconf,jconf = hconfig.gradient_config,hconfig.jacobian_config
    seeds = jconf.seeds
    duals = jconf.duals[1]
    return ForwardDiff.GradientConfig{T,V,N,typeof(duals)}(seeds,duals)
end

function _JacobianConfig(hconfig::ForwardDiff.HessianConfig{T,V,N},yduals = nothing) where {T,V,N}
    seeds = hconfig.jacobian_config.seeds
    duals = hconfig.jacobian_config.duals
    #duals = (xduals,yduals)
    #@show xduals[1]
    #@show xduals[2]
    #@show typeof(xduals)
    return ForwardDiff.JacobianConfig{T,V,N,typeof(duals)}(seeds,duals)
end

chunksize(::ForwardDiff.Chunk{C}) where {C} = C
chunksize(x::AbstractArray) = chunksize(ForwardDiff.Chunk(x))

function autochunk(x)
    k = ForwardDiff.pickchunksize(length(x))
    return ForwardDiff.Chunk{k}()
end

"""
    primalval(x::Real)

Returns the primal value of a value. Strips all duals from `ForwardDiff`. Useful in debugging:

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

@generated function primalval_struct(x::M) where M
    names = fieldnames(M)
    Base.typename(M).wrapper
    primalvals = Expr(:call,Base.typename(M).wrapper)
    for name in names
        push!(primalvals.args,:(primalval(x.$name)))
    end
    return primalvals
end

primal_eltype(x) = primal_eltype(eltype(x))
primal_eltype(::Type{W}) where W <: ForwardDiff.Dual{T,V} where {T,V} = primal_eltype(V)
primal_eltype(::Type{T}) where T = T

#eager version:
primalval_eager(x) = primalval(x)
function primalval_eager(x::AbstractArray{T}) where T <: ForwardDiff.Dual 
    res = similar(x,primal_eltype(T))
    length(res) == 0 && return res
    return map!(primalval,res,x)
end

#this struct is used to wrap a vector of ForwardDiff.Dual's and just return the primal values, without allocations
struct PrimalValVector{T,V} <: AbstractVector{T}
    vec::V
end

function PrimalValVector(v::V) where V
    T = primal_eltype(v)
    PrimalValVector{T,V}(v)
end

Base.size(x::PrimalValVector) = Base.size(x.vec)
Base.length(x::PrimalValVector) = Base.length(x.vec)
Base.@propagate_inbounds function Base.getindex(x::PrimalValVector{T},i) where T
    return primalval(x.vec[i])::T
end

#array overload for primalval
function primalval(x::AbstractArray{T}) where T <: ForwardDiff.Dual
    return PrimalValVector(x)
end
#=
gradient at index i

=#

struct GradᵢVector{T,V} <: AbstractVector{T}
    i::Int
    val::T
    vector::V
end

Base.@propagate_inbounds function Base.getindex(x::GradᵢVector{T,V},i) where {T,V}
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

sqrt_strong_zero(x) = sqrt(x)
iszeroprimal(x) = iszero(x)
iszeroprimal(x::ForwardDiff.Dual) = iszero(primalval(x))

strong_zero(y,x) = y

function strong_zero(y::ForwardDiff.Dual{T,V,P},x::ForwardDiff.Dual{T,V,P}) where {T,V,P}
    dx = ForwardDiff.partials(x)
    dy = ForwardDiff.partials(y)
    dy2 = strong_zero(dy,dx)
    return ForwardDiff.Dual{T,V,P}(x.value,dy2)
end

function strong_zero(dy::ForwardDiff.Partials{N,V},dx::ForwardDiff.Partials{N,V}) where {N,V}
    dx_tuple = dx.values
    dy_tuple = dy.values
    #is_zero_x = map(iszeroprimal,dx_tuple)
    dy2_tuple = ntuple(i -> iszeroprimal(dx_tuple[i]) ? zero(dy_tuple[i]) : dy_tuple[i],Val(N))
    return ForwardDiff.Partials{N,V}(dy2_tuple)
end

strong_zero(f::Base.Callable,x) = f(x)
strong_zero(f::Base.Callable,x::ForwardDiff.Dual) = strong_zero(f(x),x)
