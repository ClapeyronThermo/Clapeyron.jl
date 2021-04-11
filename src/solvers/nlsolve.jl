#nlsolve functionality
"""
    function nlsolve(f!,x0,method=LineSearch(Newton()),options=NEqOptions())

given a function f!(result,x) that returns a system of equations, 
`nlsolve!(f!,x0)` returns a vector x that satisfies result = 0 to some accuracy. 

uses NLSolvers.jl as backend and ForwardDiff as AD

"""


function nlsolve(f!,x0,method=LineSearch(Newton()),options=NEqOptions())
    #f! = (F,x) -> Obj_Sat(model, F, T, exp10(x[1]), exp10(x[2]))
    len = length(x0)
    xcache = zeros(eltype(x0),len)
    Fcache = zeros(eltype(x0),len)
    JCache = zeros(eltype(x0),len,len)
    function j!(J,x)
        #@show J
        ForwardDiff.jacobian!(J,f!,Fcache,x)
    end
    function fj!(F,J,x) 
        #@show J,F
        ForwardDiff.jacobian!(J,f!,F,x)
        F,J
    end
    
    function jv!(x)
        function JacV(Fv,v)
            ForwardDiff.jacobian!(JCache,f!,Fcache,v)
            Fv .= Jcache * v
        end
        return LinearMap(JacV,length(x))
    end
    vectorobj = NLSolvers.VectorObjective(f!,j!,fj!,jv!)
    vectorprob = NEqProblem(vectorobj)
    res = solve(vectorprob, x0,method , options)
end