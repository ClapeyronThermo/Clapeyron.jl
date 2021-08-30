"""
    f∂f(f,x)

returns f and ∂f/∂x evaluated in `x`, using `ForwardDiff.jl`, `DiffResults.jl` and `StaticArrays.jl` to calculate everything in one pass.
"""
function f∂f(f,x::T) where T
    _f(z) = f(only(z))
    x_vec =   SVector(x)
    ∂result = DiffResults.GradientResult(x_vec)  
    _∂f =  ForwardDiff.gradient!(∂result, _f,x_vec)
    fx =  DiffResults.value(_∂f)
    ∂f∂x = only(DiffResults.gradient(_∂f))
    return fx,∂f∂x
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
    ad_newton(f,x0;kwargs...)

    Performs newton root finding (via Roots.jl), but the derivatives are calculated with ForwardDiff.jl
    kwargs are passed to `Roots.find_zero`
"""
function ad_newton(f::F,x0;xatol=nothing, xrtol=nothing, maxevals = 100) where F
    function f0(_x)
        fx,∂f∂x = f∂f(f,_x)
        return fx,fx/∂f∂x
    end
    x = float(x0)
    T = typeof(x)
    atol = xatol !== nothing ? xatol : oneunit(T) * (eps(one(T)))^(4/5)
    rtol = xrtol !== nothing ? xrtol : eps(one(T))^(4/5)
    xo = T(Inf)
    for i in 1:maxevals
        fx, Δx = f0(x)
        iszero(fx) && return x

        x -= Δx

        if isapprox(x, xo, atol=atol, rtol=rtol)
            return x
        end
        xo = x
    end
    return zero(x0)/zero(x0)
end


