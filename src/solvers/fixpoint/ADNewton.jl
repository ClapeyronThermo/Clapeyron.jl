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


