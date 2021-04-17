#nlsolve functionality
"""
    function nlsolve(f!,x0,TrustRegion(Newton(), Dogleg()), options=NEqOptions())


given a function f!(result,x) that returns a system of equations, 
`nlsolve!(f!,x0)` returns a vector x that satisfies result = 0 to some accuracy. 

uses NLSolvers.jl as backend, the jacobian is calculated with ForwardDiff.jl

"""

function nlsolve(f!,x0,method=TrustRegion(Newton(), Dogleg()),options=NEqOptions())
    #f! = (F,x) -> Obj_Sat(model, F, T, exp10(x[1]), exp10(x[2]))
    len = length(x0)
    xcache = zeros(eltype(x0),len)
    Fcache = zeros(eltype(x0),len)
    JCache = zeros(eltype(x0),len,len)
    jconfig = ForwardDiff.JacobianConfig(f!,x0,x0)
    function j!(J,x)
        #@show J
        ForwardDiff.jacobian!(J,f!,Fcache,x,jconfig)
    end
    function fj!(F,J,x) 
        #@show J,F
        ForwardDiff.jacobian!(J,f!,F,x,jconfig)
        F,J
    end
    
    function jv!(x)
        function JacV(Fv,v)
            ForwardDiff.jacobian!(JCache,f!,Fcache,v,jconfig)
            Fv .= Jcache * v
        end
        return LinearMap(JacV,length(x))
    end
    vectorobj = NLSolvers.VectorObjective(f!,j!,fj!,jv!)
    vectorprob = NEqProblem(vectorobj)
    res = solve(vectorprob, x0,method , options)
    #@show res
end