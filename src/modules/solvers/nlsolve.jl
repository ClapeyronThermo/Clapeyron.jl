

is_rootsjl_method(method) = false

#nlsolve functionality
"""
    function nlsolve(f!,x0,method=TrustRegion(Newton(), NWI()), options=NEqOptions(),chunk = ForwardDiff.Chunk{2}())


Given a function `f!(result,x)` that returns a system of equations,
`nlsolve(f!,x0)` returns a `NLSolvers.ConvergenceInfo` struct that contains the results of the non-linear solving procedure.

Uses `NLSolvers.jl` as backend, the jacobian is calculated with `ForwardDiff.jl`, with the specified `chunk` size.

To obtain the underlying solution vector, use [`solution`](@ref).

To see available solvers and options, check `NLSolvers.jl`.
"""
function nlsolve(f!,x0,method = TrustRegion(Newton(), Dogleg()),options=NEqOptions(),chunk = ForwardDiff.Chunk{2}())
    if is_rootsjl_method(method)
        return roots_nlsolve(f!,x0,method,options)
    end
    vector_objective = ADVectorObjective(f!,x0,chunk)
    nl_problem = NEqProblem(vector_objective; inplace = __is_implace(x0))
    return nlsolve(nl_problem, x0,method, options)
end

function nlsolve(nl_problem::NEqProblem,x0,method = TrustRegion(Newton(), Dogleg()),options=NEqOptions())
    return NLSolvers.solve(nl_problem, x0,method, options)
end

function ADVectorObjective(f!,x0,chunk)
    FF = similar(x0)
    jconfig = ForwardDiff.JacobianConfig(f!,FF,x0,chunk)
    function j!(J,x)
        ForwardDiff.jacobian!(J,f!,FF,x,jconfig)
        J
    end
    function fj!(F,J,x)
        ForwardDiff.jacobian!(J,f!,F,x,jconfig)
        F,J
    end
    return NLSolvers.VectorObjective(f!,j!,fj!,nothing)
end

function ADVectorObjective(f!,x0::StaticArrays.SVector,chunk)
    f̄ = Base.Fix1(f!,nothing)
    f(F,x) = f!(nothing,x) #we assume that the F argument is unused in static arrays
    j(J,x) = ForwardDiff.jacobian(f̄,x)
    fj(F,J,x) = FJ_ad(f̄,x)
    return NLSolvers.VectorObjective(f,j,fj,nothing)
end

ADVectorObjective(f!,x0::StaticArrays.SVector) = ADVectorObjective(f!,x0,nothing)

#= only_fj!: NLsolve.jl legacy form:

function only_fj!(F, J, x)
    # shared calculations begin
    # ...
    # shared calculation end
    if !(F == nothing)
        # mutating calculations specific to f! goes here
    end
    if !(J == nothing)
        # mutating calculations specific to j! goes
    end
end
=#
function only_fj!(fj!::T) where T
    function _f!(F,x)
        fj!(F,nothing,x)
        F
    end

    function _fj!(F,J,x)
        fj!(F,J,x)
        F,J
    end

    function _j!(J,x)
        fj!(nothing,J,x)
        J
    end
    return NLSolvers.VectorObjective(_f!,_j!,_fj!,nothing) |> NEqProblem
    # return NLSolvers.VectorObjective(f!,j!,fj!,jv!) |> NEqProblem
end

#=
#trying to make nlsolve(f,x0,LineSearch(Newton(),HZAW())) work

function NLSolvers.upto_gradient(meritobj::NLSolvers.MeritObjective, ∇f, x)
    neq = meritobj.prob
    G = neq.R.F(∇f, x)
    F = (norm(G)^2) / 2
    return F,G
end
=#

#=

=#

struct Newton2Var end

function nlsolve2(f::FF,x::SVector{NN,TT},method::Newton2Var,options=NEqOptions(),tag = f;bounds = nothing) where {FF,NN,TT}
    function FJ(_z)
        return FJ_ad(f,_z)
    end
    Fx, Jx = FJ(x)
    z = x
    T = eltype(Fx)
    stoptol = T(options.f_abstol)
    ρF0, ρ2F0 = norm(Fx, Inf), norm(Fx, 2)
    nan = T(NaN)
    ρs = nan
    #@show ρF0
    if ρF0 < stoptol
        return x
    end
    iter = 1
    converged = false
    while iter ≤ options.maxiter
        d = Jx \ -Fx
        #@show Jx, Fx
        x = __new_x_newton2var(x,d,bounds)
        Fx, Jx = FJ(x)
        in_bounds = __newton2var_in_bounds(x,d,bounds)
        ρF = bounded_norm(Fx, Inf, in_bounds)
        ρs = bounded_norm(Fx, Inf, in_bounds)
        #@show ρF, ρs
        if ρs <= stoptol || ρF <= stoptol
            converged = true
            break
        end

        if !all(isfinite,x)
            converged = false
            break
        end
        iter += 1
    end
    if !converged
        x  = nan .* x
    end
    return x
end



bounded_norm(x,nn,in_bounds) = norm(x .* in_bounds,nn)
bounded_norm(x,nn,::Nothing) = norm(x,nn)

__newton2var_in_bounds(x,d,::Nothing) = nothing
function __newton2var_in_bounds(x,d,bounds)
    lb,ub = bounds
    return (lb .<= (x + d) .<= ub) 
end

function __new_x_newton2var(x::A,d::B,::Nothing) where {A,B}
    return x + d
end

function __new_x_newton2var(x::A,d::B,bounds) where {A,B}
    lb,ub = bounds
    return bound_1d.(x,d,lb,ub)
end

function bound_1d(x,d,lb,ub)
    y = x + d
    if y < lb
        return 0.5*(x + lb)
    elseif y > ub
        return 0.5*(x + ub)
    else
        return y
    end
end
#=
Roots.jl extension
=#

is_rootsjl_method(method::Roots.AbstractUnivariateZeroMethod) = true

function roots_nlsolve(f::F,x0,method::Roots.AbstractBracketingMethod,options) where F
    prob = Roots.ZeroProblem(f,x0)
    brk = f(x0[1])*f(x0[2])
    if brk > 0
        sol = zero(brk)/zero(brk)
    else
        sol = Roots.solve(prob,method)
    end
end

function roots_nlsolve(f::F,x0::Number,method::Roots.AbstractNonBracketingMethod ,options) where F
    prob = Roots.ZeroProblem(f,x0)
    sol = Roots.solve(prob,method)
end

function roots_nlsolve(f::F,x0::Number,method::Roots.AbstractNewtonLikeMethod ,options) where F
    prob = Roots.ZeroProblem(to_newton(f),x0)
    sol = Roots.solve(prob,method)
end

function roots_nlsolve(f::F,x0::Number,method::Roots.AbstractHalleyLikeMethod,options) where F
    prob = Roots.ZeroProblem(to_halley(f),x0)
    sol = Roots.solve(prob,method)
end

#iterative solver

function solve1_update_state(state, x, fx; full_iter=true)
    xs, fxs, damp, status = state
    x1, x2, x3 = xs
    f1, f2, f3 = fxs

    if status == :no_init
        return ((x, x2, x3), (fx, f2, f3), damp, :iter0)

    elseif status == :iter0
        return ((x, x1, x2), (fx, f1, f2), damp, :iter_initial)

    elseif status == :iter_initial
        new_damp = max(zero(damp), damp - 0.2)
        new_status = iszero(new_damp) ? :iter_full : :iter_initial
        return ((x, x1, x2), (fx, f1, f2), new_damp, new_status)

    elseif status == :iter_full
        # Only attempt to establish a bracket when fx is a reliable evaluation
        if full_iter
            if fx * f1 < 0
                na, nb, nfa, nfb = Roots.sort_smallest(x, x1, fx, f1)
                # c: best of the two remaining points
                nc, nfc = abs(f2) <= abs(f3) ? (x2, f2) : (x3, f3)
                return ((na, nb, nc), (nfa, nfb, nfc), damp, :bounded)
            elseif fx * f2 < 0
                na, nb, nfa, nfb = Roots.sort_smallest(x, x2, fx, f2)
                nc, nfc = abs(f1) <= abs(f3) ? (x1, f1) : (x3, f3)
                return ((na, nb, nc), (nfa, nfb, nfc), damp, :bounded)
            elseif fx * f3 < 0
                na, nb, nfa, nfb = Roots.sort_smallest(x, x3, fx, f3)
                nc, nfc = abs(f1) <= abs(f2) ? (x1, f1) : (x2, f2)
                return ((na, nb, nc), (nfa, nfb, nfc), damp, :bounded)
            end
        end

        # No bracket found, or full_iter=false: maintain 3 best points sorted by x.
        _, imax = findmax((abs(f1),abs(f2),abs(f3)))
        keep = if imax == 1
            ((x2, f2), (x3, f3), (x, fx))
        elseif imax == 2
            ((x1, f1), (x3, f3), (x, fx))
        else
            ((x1, f1), (x2, f2), (x, fx))
        end
        w1, w2, w3 = sort(keep)
        return ((w1[1], w2[1], w3[1]), (w1[2], w2[2], w3[2]), damp, :iter_full)

    elseif status == :bounded
        a, b, c = x1, x2, x3
        fa, fb, fc = f1, f2, f3
        # Invariant: fa*fb < 0, |fa| <= |fb|, c is the previously displaced endpoint.

        if !full_iter
            # fx is approximate: don't risk corrupting the bracket.
            # Only update c if the new point is strictly better, otherwise hold.
            if abs(fx) < abs(fc)
                return ((a, b, x), (fa, fb, fx), damp, :bounded)
            else
                return state
            end
        end

        # Fix: correctly track which endpoint gets displaced.
        # The displaced endpoint (the one leaving the bracket) becomes the new c,
        # which is what gives Chandrapatla's ξ and ϕ their meaning.
        if fx * fa < 0
            # New bracket: (x, a). Old b leaves → new c = b.
            na, nb, nfa, nfb = Roots.sort_smallest(x, a, fx, fa)
            return ((na, nb, b), (nfa, nfb, fb), damp, :bounded)
        elseif fx * fb < 0
            # New bracket: (x, b). Old a leaves → new c = a.
            na, nb, nfa, nfb = Roots.sort_smallest(x, b, fx, fb)
            return ((na, nb, a), (nfa, nfb, fa), damp, :bounded)
        else
            # fa*fb < 0 by invariant, so x must bracket with one side.
            # Reaching here means fx ≈ 0 or a noisy evaluation slipped through.
            # Guard: update c only if the new point is better than what we have.
            if abs(fx) < abs(fc)
                return ((a, b, x), (fa, fb, fx), damp, :bounded)
            else
                return state
            end
        end

    else
        return state
    end
end


function solve1_new_iter(old_state,x,fx,dfx = nothing;full_iter = true)
    ∂fx = dfx == nothing ? fx : oftype(fx,dfx)
    state = solve1_update_state(old_state,x,fx;full_iter)
    xx,fxx,α,status = state
    if status == :iter0
        if dfx == nothing
            x0 = x
        else
            x0 = (x - (1 - α)*f/∂fx)
        end
        return x0,state
    elseif status == :iter_initial
        xa,xb,_ = xx
        fa,fb,_ = fxx
        dfdx = dfx == nothing ? ((fb - fa)/(xb - xa)) : ∂fx   
        xnew = (xa - (1 - α)*fa/dfdx)
        return xnew,state
    elseif status == :iter_full
        x1, x2, x3 = xx
        f1, f2, f3 = fxx
        if full_iter
            xnew = Roots.inverse_quadratic_step(x1, x2, x3, f1, f2, f3)
            xlo, xhi = extrema((x1,x2,x3))
            span = xhi - xlo
            if !(xlo - span < xnew < xhi + span)
                dfdx = dfx == nothing ? ((f2 - f1) / (x2 - x1)) : ∂fx
                xnew = x1 - f1 / dfdx
            end
        else
            dfdx = dfx == nothing ? ((f2 - f1) / (x2 - x1)) : ∂fx
            xnew = x1 - f1 / dfdx
        end

        if f1*f2 < 0
            if x1 < x2
                xlo,xhi,flo,fhi = x1,x2,f1,f2
            else
                xlo,xhi,flo,fhi = x2,x1,f2,f1
            end
            if !(xlo <= xnew <= xhi)
                xm = 0.5*(xlo + xhi)
                xrf = (fb * a - fa * b) / (fb - fa)

            if xlo < xrf < xhi
                # Check Illinois adjustment: if regula falsi keeps one side fixed,
                # the Illinois half-weight step can do better.
                fa_ill = fa / 2
                xill = ((fb) * a - fa_ill * b) / (fb - fa_ill)
                xnew = (xlo < xill < xhi) ? xill : xrf
            else
                xnew = xm  # regula falsi left the bracket, bisect
            end
            end
        end

        return xnew, state
    elseif status == :bounded
        a, b, c = xx
        fa, fb, fc = fxx
        # a is best (|fa| <= |fb|), [a,b] is the bracket, c is displaced endpoint.

        xlo, xhi = minmax(a, b)
        xm = (a + b) / 2  # bisection

        # Chandrapatla's condition: IQI is interpolating (not extrapolating).
        ξ = (a - b) / (c - b)
        ϕ = (fa - fb) / (fc - fb)
        ϕ² = ϕ^2
        Δ = (ϕ² < ξ) && (1 - 2ϕ + ϕ² < 1 - ξ)

        if Δ
            xnew = Roots.inverse_quadratic_step(a, b, c, fa, fb, fc)
            # Guard: IQI must land strictly inside the bracket.
            if !(xlo < xnew < xhi)
                xnew = xm
            end
        else
            # Illinois step: modify the regula falsi weight at the retained endpoint
            # to avoid one-sided stagnation, then guard it. Falls back to bisection.
            #
            # Standard regula falsi: xrf = (fb*a - fa*b) / (fb - fa)
            # Illinois: halve the function value at whichever endpoint is retained.
            # If the new point lands on the same side as a (the best point),
            # the retained endpoint is b, so we use fb/2 in the formula.
            xrf = (fb * a - fa * b) / (fb - fa)

            if xlo < xrf < xhi
                # Check Illinois adjustment: if regula falsi keeps one side fixed,
                # the Illinois half-weight step can do better.
                fa_ill = fa / 2
                xill = ((fb) * a - fa_ill * b) / (fb - fa_ill)
                xnew = (xlo < xill < xhi) ? xill : xrf
            else
                xnew = xm  # regula falsi left the bracket, bisect
            end
        end
        return xnew, state
    else
        return first(xx)*NaN,state
    end
    return first(xx)*NaN,state
end