

function ADScalarObjective(f,x0::T) where T <: AbstractArray
    Hresult = DiffResults.HessianResult(x0)
    Hconfig = ForwardDiff.HessianConfig(f,Hresult,x0)
    function fg(_df,x)
        Gconfig = _GradientConfig(Hconfig)
        result = ForwardDiff.gradient!(Hresult,f,x,Gconfig)
        if __is_implace(x)
            _df .= DiffResults.gradient(result)
            df = _df
        else
            df = DiffResults.gradient(result)
        end
        fx = DiffResults.value(result)
        return fx,df
    end

      function g(df,x)
        _,df = fg(df,x)
        return df
    end

    function fgh(df,d2f,x)
        result = ForwardDiff.hessian!(Hresult,f,x,Hconfig)
        fx = DiffResults.value(result)
        if __is_implace(x)
            d2f .= DiffResults.hessian(result)
            df .= DiffResults.gradient(result)
            return fx,df,d2f
        else
            _d2f = DiffResults.hessian(result)
            _df = DiffResults.gradient(result)
            return fx,_df,_d2f
        end 
    end

    function h(d2f,x)
        result = ForwardDiff.hessian!(Hresult,f,x)
        if __is_implace(x)
            d2f .= DiffResults.hessian(result)
            return d2f
        else
            return DiffResults.hessian(result)
        end
    end

    return ScalarObjective(f=f,
    g=g,
    fg=fg,
    fgh=fgh,
    h=h)
end

function ADScalarObjective(f::F,x0::T) where {F,T <: SVector}
    function g(_a1,x)
        return ForwardDiff.gradient(f,x)
    end
    
    function fg(_a1,x)
        gres = DiffResults.GradientResult(x)
        result = ForwardDiff.gradient!(gres,f,x)
        df = DiffResults.gradient(result)
        fx = DiffResults.value(result)
        return fx,df
    end

    function fgh(_a1,_a2,x)
        Hresult = DiffResults.HessianResult(x)
        result = static_fgh(Hresult,f,x)
        fx = DiffResults.value(result)
        _d2f = DiffResults.hessian(result)
        _df = DiffResults.gradient(result)
        return fx,_df,_d2f
    end

    function h(_a1,x)
        return ForwardDiff.hessian(f,x)
    end

    return ScalarObjective(f=f,
    g=g,
    fg=fg,
    fgh=fgh,
    h=h)
end


function ADScalarObjective(f,x0::Number)
    function g(x)
        return derivative(f,x)
    end

    function fg(_,x)
        ∂fx,x = f∂f(f,x)
        return ∂fx,x
    end

    function fgh(_,_,x)
        fx,∂fx,∂2fx = f∂f∂2f(f,x)
        return fx,∂fx,∂2fx
    end

    function h(_,x)
        ∂2fx = derivative(g,x)
        return ∂2fx
    end

    return ScalarObjective(f = f,
    g = g,
    fg = fg,
    fgh = fgh,
    h = h)
end

#uses brent, the same default that Optim.jl uses
function optimize(f,x0::NTuple{2,T},method=BrentMin(),options=OptimizationOptions()) where T<:Real
    scalarobj = ADScalarObjective(f,first(x0))
    optprob = OptimizationProblem(obj = scalarobj,bounds = x0, inplace=false)
    res = NLSolvers.solve(optprob,method,options)
    return res
end

#general one, with support for ActiveBox
function optimize(f,x0,method=LineSearch(Newton2(x0),NLSolvers.Static(1.0)),options=OptimizationOptions();bounds = nothing)
    inplace = __is_implace(x0)
    scalarobj = ADScalarObjective(f,x0)
    optprob = OptimizationProblem(obj = scalarobj,inplace = inplace,bounds = bounds)
    return NLSolvers.solve(optprob,x0,method,options)
end

function optimize(optprob::OptimizationProblem,x0,method=LineSearch(Newton2(x0),NLSolvers.Static(1.0)),options=OptimizationOptions();bounds = nothing)
    return NLSolvers.solve(optprob,x0,method,options)
end

#build scalar objective -> Optimization Problem
function optimize(scalarobj::ScalarObjective,x0,method=LineSearch(Newton2(x0),NLSolvers.Static(1.0)),options=OptimizationOptions();bounds = nothing)
    inplace = __is_implace(x0)
    optprob = OptimizationProblem(obj = scalarobj,inplace = inplace,bounds = bounds)
    return NLSolvers.solve(optprob,x0,method,options)
end

function optimize(f,x0,method::NLSolvers.NelderMead,options=OptimizationOptions();bounds = nothing)
    scalarobj = ScalarObjective(f = f)
    optprob = OptimizationProblem(obj = scalarobj,inplace = false,bounds = bounds)
    return NLSolvers.solve(optprob,x0,method,options)
end


#= only_fg!: Optim.jl legacy form:
function fg!(F,G,x)
  # do common computations here
  # ...
  if G != nothing
    # code to compute gradient here
    # writing the result to the vector G
  end
  if F != nothing
    # value = ... code to compute objective function
    return value
  end
end
=#

function only_fg!(fg!::T) where T
    function f(x)
        return fg!(true,nothing,x)
    end

    function g(df,x)
        fg!(nothing,df,x)
        return df
    end
    function fg(df,x)
        fx = fg!(true,df,x)
        return fx,df
    end

    return ScalarObjective(f=f,
    g=g,
    fg=fg,
    fgh=nothing,
    h=nothing)
end

#= only_fgh!: Optim.jl legacy form:
function fgh!(F,G,H,x)
  G == nothing || # compute gradient and store in G
  H == nothing || # compute Hessian and store in H
  F == nothing || return f(x)
  nothing
end
=#

function only_fgh!(fgh!::T) where T
    function f(x)
        return fgh!(true,nothing,nothing,x)
    end

    function g(df,x)
        fgh!(nothing,df,nothing,x)
        return df
    end

    function fg(df,x)
        fx = fgh!(true,df,nothing,x)
        return fx,df
    end

    function fgh(df,d2f,x)
        fx = fgh!(true,df,d2f,x)
        return fx,df,d2f
    end

    function h(d2f,x)
        fgh!(nothing,nothing,d2f,x)
        return d2f
    end

    return ScalarObjective(f=f,
    g=g,
    fg=fg,
    fgh=fgh,
    h=h)
end

struct BoundOptim1Var end

function optimize(f,x0::NTuple{2,T},method::BoundOptim1Var,options=OptimizationOptions()) where T<:Real
    return _1var_optimize_quad(f,x0)
end

quad_interp(x,f) = quad_interp(x[1],x[2],x[3],f[1],f[2],f[3])

function quad_interp(xa,xb,xc,fa,fb,fc)
        #f1 = ax12 + bx1 + c
    #f2 = ax22 + bx2 + c
    #f3 = ax32 + bx3 + c
    A = @SMatrix [xa*xa xa oneunit(xa); xb*xb xb oneunit(xb); xc*xc xc oneunit(xc)]
    B = SVector((fa,fb,fc))
    z = A\B
    a,b,c = z
    return a,b,c
end

function _1var_optimize_quad(f,x0)
    xa,xb = minmax(x0[1],x0[2])
    xa0,xb0 = xa,xb
    fa,fb = f(xa),f(xb)
    xc = 0.5*(xa + xb)
    fc = f(xc)
    for i in 1:20
        a,b,c = quad_interp(xa,xb,xc,fa,fb,fc)
        xmin = -b/(2*a)
        fmin = f(xmin)
        _f = SVector((fa,fb,fc,fmin))
        f_worst,idx = findmax(_f)
        if fmin != f_worst && isfinite(fmin)
            _x = SVector((xa,xb,xc,xmin))
            xa,xb,xc = StaticArrays.deleteat(_x,idx)
            fa,fb,fc = StaticArrays.deleteat(_f,idx)
        elseif fmin == f_worst
            _f2 = SVector((fa,fb,fc))
            _x2 = SVector((xa,xb,xc))
            _,idx2 = findmax(_f2)
            xa,xb = StaticArrays.deleteat(_x2,idx2)
            fa,fb = StaticArrays.deleteat(_f2,idx2)
            xc = 0.5*(xa + xb)
            fc = f(xc)
        else
            return zero(fmin)/zero(fmin)
        end
        xmin,xmax = extrema((xa,xb,xc))
        fxx = minimum((fa,fb,fc))
        #@show abs(xmin - xmax),fxx
        if abs(xmin - xmax) < sqrt(eps(xmin))
            break
        end
    end
    xmin,xmax = extrema((xa,xb,xc))
    #@show xa0,xmin,xmax,xb0
    return 0.5*(xmin + xmax)
end