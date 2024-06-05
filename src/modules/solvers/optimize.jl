
function ADScalarObjective(f,x0::AbstractArray,chunk = autochunk(x0),val::Val{N} = Val{2}()) where N
    Hres = DiffResults.HessianResult(x0)
    function _g(df,x,Gresult)
        ForwardDiff.gradient!(Gresult,f,x)
        df .= DiffResults.gradient(Gresult)
        df
    end
    
    function _fg(df,x,Gresult)
        ForwardDiff.gradient!(Gresult,f,x)
        df .= DiffResults.gradient(Gresult)
        fx = DiffResults.value(Gresult)
        return fx,df
    end

    function _fgh(df,d2f,x,Hresult)
        ForwardDiff.hessian!(Hresult,f,x)
        d2f .= DiffResults.hessian(Hresult)
        df .= DiffResults.gradient(Hresult)
        fx = DiffResults.value(Hresult)
        return fx,df,d2f
    end

    function h(d2f,x)
        ForwardDiff.hessian!(d2f,f,x)
        d2f
    end
    g(df,x) = _g(df,x,Hres)
    fg(df,x) = _fg(df,x,Hres)
    fgh(df,d2f,x) = _fgh(df,d2f,x,Hres)
    return ScalarObjective(f=f,
    g=g,
    fg=fg,
    fgh=fgh,
    h=h)
end


function ADScalarObjective(f,x0::Number,autochunk)
    function g(x)
        return derivative(f,x)
    end

    function fg(∂fx,x)
        ∂fx,x = f∂f(f,x)
        return ∂fx,x
    end

    function fgh(∂fx,∂2fx,x)
        fx,∂fx,∂2fx = f∂f∂2f(f,x)
        return fx,∂fx,∂2fx
    end

    function h(∂2fx,x)
        ∂2fx = derivative(g,x)
        return ∂2fx
    end

    return ScalarObjective(f=f,
    g=g,
    fg=fg,
    fgh=fgh,
    h=h)
end
#uses brent, the same default that Optim.jl uses
function optimize(f,x0::NTuple{2,T},method=BrentMin(),options=OptimizationOptions()) where T<:Real
    scalarobj = ADScalarObjective(f,first(x0),nothing)
    optprob = OptimizationProblem(obj = scalarobj,bounds = x0, inplace=false)
    res = NLSolvers.solve(optprob,method,options)
    return res
end
#general one, with support for ActiveBox
function optimize(f,x0,method=LineSearch(Newton()),options=OptimizationOptions();bounds = nothing)
    scalarobj = ADScalarObjective(f,x0,autochunk)
    optprob = OptimizationProblem(obj = scalarobj,inplace = (x0 isa Number),bounds = bounds)
    return NLSolvers.solve(optprob,x0,method,options)
end

function optimize(optprob::OptimizationProblem,x0,method=LineSearch(Newton()),options=OptimizationOptions();bounds = nothing)
    return NLSolvers.solve(optprob,x0,method,options)
end
#build scalar objective -> Optimization Problem
function optimize(scalarobj::ScalarObjective,x0,method=LineSearch(Newton()),options=OptimizationOptions();bounds = nothing)
    optprob = OptimizationProblem(obj = scalarobj,inplace = (x0 isa Number),bounds = bounds)
    return NLSolvers.solve(optprob,x0,method,options)
end

function optimize(f,x0,method::NLSolvers.NelderMead,options=OptimizationOptions();bounds = nothing)
    scalarobj = ScalarObjective(f = f)
    optprob = OptimizationProblem(obj = scalarobj,inplace = (x0 isa Number),bounds = bounds)
    return NLSolvers.solve(optprob,x0,method,options)
end

x_minimum(res::NLSolvers.ConvergenceInfo) = res.info.minimum
#for BrentMin (should be fixed at NLSolvers 0.3)
x_minimum(res::Tuple{<:Number,<:Number}) = last(res)

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
