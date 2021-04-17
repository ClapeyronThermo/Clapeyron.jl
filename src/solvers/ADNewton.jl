"""
    autonewton(f,x)

returns f/(df/fx) evaluated in `x`, using `ForwardDiff.jl`, `DiffResults.jl` and `StaticArrays.jl` to calculate f and dfdx in one pass
"""
function autonewton(f,x::T) where T
    _f(z) = f(only(z))
    x_vec =   SVector(x)
    ∂result = DiffResults.GradientResult(x_vec)  
    _∂f =  ForwardDiff.gradient!(∂result, _f,x_vec)
    fx =  DiffResults.value(_∂f)
    ∂f∂x = only(DiffResults.gradient(_∂f))
    return fx,fx/∂f∂x
end

"""
    ad_newton(f,x0;kwargs...)

    Performs newton root finding (via Roots.jl), but the derivatives are calculated with ForwardDiff.jl
    kwargs are passed to `Roots.find_zero`
"""
function ad_newton(f,x0;kwargs...)
    f0 = x-> autonewton(f,x)
    return Roots.find_zero(f0,x0,Roots.Newton();kwargs...)
end

