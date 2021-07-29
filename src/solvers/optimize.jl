"""
    function optimize(f,x0,method=LineSearch(Newton()), options=OptimizationOptions())


"""
function optimize(f,x0::Number,method=LineSearch(Newton()),options=OptimizationOptions())

    function g(x)
        return ForwardDiff.derivative(f,x)
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
        ∂2fx = ForwardDiff.derivative(g,x)
        return ∂2fx
    end

    scalarobj = ScalarObjective(f=f,
    g=g,
    fg=fg,
    fgh=fgh,
    h=h)
    
    optprob = OptimizationProblem(scalarobj; inplace=false) 
    return NLSolvers.solve(optprob, x0, method,options)
end