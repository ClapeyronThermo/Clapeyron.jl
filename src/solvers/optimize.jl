
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

#=
function ADScalarObjective(f,x0::Number)
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
"""
    function optimize(f,x0,method=LineSearch(Newton()), options=OptimizationOptions())
"""
=#

function optimize(f,x0,method=LineSearch(Newton()),chunk =autochunk(x0),options=OptimizationOptions())
    scalarobj = ADScalarObjective(f,x0,chunk)   
    optprob = OptimizationProblem(scalarobj; inplace=false) 
    return NLSolvers.solve(optprob, x0, method,options)
end

x_minimum(res::NLSolvers.ConvergenceInfo) = res.info.minimum