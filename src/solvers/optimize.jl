
function ADScalarObjective(f,x0::AbstractArray,chunk = autochunk(x0))
    Hres = DiffResults.HessianResult(x0)
    function _g(df,x,Hresult)
        ForwardDiff.gradient!(Hresult,f,x)
        df .= DiffResults.gradient(Hresult)
        df
    end
    function _fg(df,x,Hresult)
        ForwardDiff.gradient!(Hresult,f,x)
        df .= DiffResults.gradient(Hresult)
        fx = DiffResults.value(Hresult)
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
    fgh(df,d2f,x) =  _fgh(df,d2f,x,Hres)
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
function optimize(f,x0::NTuple{T,T},method=BrentMin(T((3 - sqrt(5)) / 2)),options=OptimizationOptions()) where T<:Real
    scalarobj = ADScalarObjective(f,x0)   
    optprob = OptimizationProblem(scalarobj;bounds = x0, inplace=false) 
    return NLSolvers.solve(optprob,method,options)
end
#general one, with support for ActiveBox
function optimize(f,x0,method=LineSearch(Newton()),options=OptimizationOptions();bounds = nothing)
    scalarobj = ADScalarObjective(f,x0,autochunk)   
    optprob = OptimizationProblem(scalarobj,inplace = (x0 isa number),bounds = bounds) 
    return NLSolvers.solve(optprob,x0,method,options)
end

function optimize(optprob::OptimizationProblem,method=LineSearch(Newton()),options=OptimizationOptions();bounds = nothing)
    return NLSolvers.solve(optprob,x0,method,options)
end

x_minimum(res::NLSolvers.ConvergenceInfo) = res.info.minimum