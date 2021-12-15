"""
    f∂f(f,x)

returns f and ∂f/∂x evaluated in `x`, using `ForwardDiff.jl`, `DiffResults.jl` and `StaticArrays.jl` to calculate everything in one pass.
"""
function f∂f(f::F, x::R) where {F,R}
    T = typeof(ForwardDiff.Tag(f, R))
    out = f(ForwardDiff.Dual{T}(x, one(x)))
    return ForwardDiff.value(out), ForwardDiff.extract_derivative(T, out)
end

"""
    f∂f∂2f(f,x)

returns f,∂f/∂x,and ∂²f/∂²x and evaluated in `x`, using `ForwardDiff.jl`, `DiffResults.jl` and `StaticArrays.jl` to calculate everything in one pass.
"""
function f∂f∂2f(f::F,x::T) where {F,T}
    _f(z) = f(only(z))
    x_vec =   SVector(x)
    ∂result = DiffResults.HessianResult(x_vec)  
    _∂f =  ForwardDiff.hessian!(∂result, _f,x_vec)
    fx =  DiffResults.value(_∂f)
    ∂f∂x = only(DiffResults.gradient(_∂f))
    ∂²f∂²x =  only(DiffResults.hessian(_∂f))
    return fx,∂f∂x,∂²f∂²x
end



function newton_fixpoint(f::F,x::T,atol) where {F,T}
    fx, gx = f(x)
    abs(fx) < atol && return x
    Δx = fx/gx
    x =x - Δx
    return x
end

"""
    ad_newton(f,x0;atol = oneunit(T) * (eps(one(T)))^(4/5),rtol = )

    Performs newton root finding (via Roots.jl), but the derivatives are calculated with ForwardDiff.jl
    kwargs are passed to `Roots.find_zero`
"""
function ad_newton(f::F,x0::T;
    rtol= (eps(one(T)))^(4/5),
    atol=rtol*oneunit(T),
    max_iters = 100) where {F,T}
    fdf(x) =  f∂f(f,x)
    f0(x) = newton_fixpoint(fdf,x,atol)
    return fixpoint(f0,x0;rtol,atol,max_iters)
end

function newton(f::F,x0::T;
    rtol= (eps(one(T)))^(4/5),
    atol=rtol*oneunit(T),
    max_iters = 100) where {F,T}
    f0(x) = newton_fixpoint(f,x,atol)
    return fixpoint(f0,x0;rtol,atol,max_iters)
end

function halley_fixpoint(f::F,x::T,atol) where {F,T}
    ff,gg,hh = f(x)
    abs(ff) < atol && return x
    Δx = ff/gg/(1-ff*hh/(2*gg^2))
    x =x - Δx
    return x
end

function halley(f::F,x0::T;
    rtol= (eps(one(T)))^(4/5),
    atol=rtol*oneunit(T),
    max_iters = 100) where {F,T}
    f0(x) = halley_fixpoint(f,x,atol)
    return fixpoint(f0,x0;rtol,atol,max_iters)
end


