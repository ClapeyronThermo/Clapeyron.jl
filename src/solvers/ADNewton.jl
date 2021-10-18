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
function f∂f∂2f(f,x::T) where T
    _f(z) = f(only(z))
    x_vec =   SVector(x)
    ∂result = DiffResults.HessianResult(x_vec)  
    _∂f =  ForwardDiff.hessian!(∂result, _f,x_vec)
    fx =  DiffResults.value(_∂f)
    ∂f∂x = only(DiffResults.gradient(_∂f))
    ∂²f∂²x =  only(DiffResults.hessian(_∂f))
    return fx,∂f∂x,∂²f∂²x
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
    return newton(Base.Fix1(f∂f,f),x0;rtol,atol,max_iters)
end

function newton(f::F,x0::T;
    rtol= (eps(one(T)))^(4/5),
    atol=rtol*oneunit(T),
    max_iters = 100) where {F,T} 
    xo = T(Inf)
    x = x0
    for _ in 1:max_iters
        fx, gx = f(x)
        Δx = fx/gx
        abs(fx) < eps(x) && return x
        x -= Δx
        if isapprox(x, xo, atol=atol, rtol=rtol)
            return x
        end
        xo = x
    end
    return zero(x0)/zero(x0)
end

function halley(fgh::F,x0::T;
    rtol= (eps(one(T)))^(4/5),
    atol=rtol*oneunit(T),
    max_iters = 100) where {F,T}
    for _ in 1:max_iters
        ff,gg,hh = fgh(x0)
        abs(ff) < eps(x0) && return x0      
        d = ff/gg/(1-ff*hh/(2*gg^2))
        if isapprox(x0, x0-d, atol=atol, rtol=rtol)
            return x0
        end
        x0=x0-d
    end
    return zero(x0)/zero(x0)
end

