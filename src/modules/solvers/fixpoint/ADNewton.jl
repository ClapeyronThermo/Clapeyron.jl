function newton_fixpoint(f::F,x::T,atol) where {F,T}
    fx, gx = f(x)
    abs(fx) < atol && return x
    Δx = fx/gx
    x =x - Δx
    return x
end

"""
    ad_newton(f,x0;atol = oneunit(T) * (eps(one(T)))^(4/5),rtol = )

    Performs Newton's root finding (via Roots.jl), but the derivatives are calculated with ForwardDiff.jl,
    kwargs are passed to `Roots.find_zero`.
"""
function ad_newton(f::F,x0::T;
    rtol= (eps(one(T)))^(4/5),
    atol=rtol*oneunit(T),
    max_iters = 100) where {F,T}
    fdf = f∂f(f)
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

function to_newton(f,x)
    f,df = f∂f(f,x)
    return f,f/df
end

to_newton(f::F) where F = Base.Fix1(to_newton,f)

function to_halley(f,x)
    f,df,d2f = f∂f∂2f(f,x)
    return f,f/df,df/d2f
end

to_halley(f::F) where F = Base.Fix1(to_halley,f)
