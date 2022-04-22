#nlsolve functionality
"""
    function nlsolve(f!,x0,method=TrustRegion(Newton(), NWI()),chunk = ForwardDiff.Chunk{2}() options=NEqOptions())


Given a function `f!(result,x)` that returns a system of equations,
`nlsolve(f!,x0)` returns a `NLSolvers.ConvergenceInfo` struct that contains the results of the non-linear solving procedure.

Uses `NLSolvers.jl` as backend, the jacobian is calculated with `ForwardDiff.jl`, with the specified `chunk` size

To obtain the underlying solution vector, use [`x_sol`](@ref)

To see available solvers and options, check `NLSolvers.jl`
"""
function nlsolve(f!,x0,method=TrustRegion(Newton(), NWI()),chunk = ForwardDiff.Chunk{2}(),options=NEqOptions();)
    vectorobj = autoVectorObjective(f!,x0,chunk)
    vectorprob = NEqProblem(vectorobj)
    return NLSolvers.solve(vectorprob, x0,method, options)
end

function autoVectorObjective(f!,x0,chunk)
    Fcache = x0 .* false
    jconfig = ForwardDiff.JacobianConfig(f!,x0,x0,chunk)
    function j!(J,x)
        ForwardDiff.jacobian!(J,f!,Fcache,x,jconfig)
    end
    function fj!(F,J,x)
        ForwardDiff.jacobian!(J,f!,F,x,jconfig)
        F,J
    end
    function jv!(x)
        function JacV(dy,v)
            return jacvec!(dy,f!,x,v)
        end
        return LinearMap(JacV,length(x))
    end
    return NLSolvers.VectorObjective(f!,j!,fj!,jv!)
end

#from SparseDiffTools.jl, but it happens to work on dense vectors as well

struct DeivVecTag end

function jacvec!(dy, f, x, v,
                      cache1 = ForwardDiff.Dual{DeivVecTag}.(x, v),
                      cache2 = ForwardDiff.Dual{DeivVecTag}.(x, v))
    f(cache2,cache1)
    dy .= ForwardDiff.partials.(cache2, 1)
end

function jacvec(f, x, v)
    partials.(f(Dual{DeivVecTag}.(x, v)), 1)
end
