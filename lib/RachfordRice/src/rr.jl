"""
    RachfordRiceResult{T}

Result returned by [`rachfordrice`](@ref).

# Fields
- `V::T`: molar vapour fraction `β`
- `L::T`: molar liquid fraction `1 - β`
- `status::RachfordRiceStatus`: phase state found, one of `RRLiquid`, `RRVapour`, `RREq`, `RRTrivial`, `RRFailure` (see [`RachfordRiceStatus`](@ref)).
- `negative_flash::Bool`: whether the problem was solved allowing `β` outside `[0, 1]`.
- `solved::Bool`: `true` if a root of the Rachford-Rice equation was found by the solver (or recovered by the boundary projection). It is `false` for single-phase shortcuts, trivial solutions and failures.
- `feasible::Bool`: `true` if the K-values straddle one (`minimum(K) <= 1 <= maximum(K)`) and have a valid domain, i.e. a two-phase root can exist.
- `iters::Int`: number of iterations used by the iterative solver (`0` when the closed-form solution, a shortcut was used or no solution is possible).

When no valid result is available, `V` and `L` are `NaN`. For a single-phase result with `negative_flash = false`, `(V, L)` is `(0, 1)` for `RRLiquid` or `(1, 0)` for `RRVapour`, with `solved = false`.
"""
struct RachfordRiceResult{T}
    V::T
    L::T
    status::RachfordRiceStatus
    negative_flash::Bool
    solved::Bool
    feasible::Bool
    iters::Int
end

function RachfordRiceResult{T}() where T
    nan = T(NaN)
    return RachfordRiceResult{T}(nan,nan,RRFailure,false,false,false,0)
end


"""
    rachfordrice(K, z, β0 = nothing; negative_flash = false, K_tol = 4*eps(eltype(K)),
                 non_inx = nothing, non_iny = nothing, verbose = nothing)
    rachfordrice(RR::RachfordRiceResidual, β0 = nothing; negative_flash = false,
                 K_tol = 4*eps(eltype(RR)), verbose = nothing)

Solve the Rachford-Rice mass balance for the molar vapour fraction `β`:

    Σᵢ zᵢ (Kᵢ - 1) / (1 + β (Kᵢ - 1)) = 0

and return a [`RachfordRiceResult`](@ref).

The solver first detects trivial, failed and single-phase cases, then computes the Whitson-Michelsen bounds on `β`. 
Small problems (up to 3 non-trivial components) are solved in closed form; otherwise a bracketed Halley iteration is used. 
When `negative_flash = false` and the root lies within `K_tol^(1/3)` of `0` or `1`, a boundary check may project the result onto the bubble point (`β = 0`) or the dew point (`β = 1`).


- `K`: vector of K-values `Kᵢ = yᵢ/xᵢ`. They must be non-negative; `Inf` and `0` are allowed.
- `z`: vector of overall molar composition. It does not need to be normalised.
- `β0`: optional initial guess for the vapour fraction. It is ignored if it falls outside the valid bounds.

`K`, `z`, `non_inx`, `non_iny` are arguments used to instantiate a [`RachfordRiceResidual`](@ref).

# Keyword arguments
- `negative_flash`: if `true`, `β` may be outside `[0, 1]` (negative flash). If `false` (default), single-phase cases return `β = 0` (liquid) or `β = 1` (vapour) and the result is clamped to `[0, 1]`.
- `K_tol`: tolerance for the trivial-solution check (`norm(K,Inf) <= K_tol`) and for the boundary tests.
- `verbose`: `nothing` (default) or a `Bool`. If `true`, diagnostic information is logged with `@info`.

# Returns
A [`RachfordRiceResult`](@ref).

# Example
```julia
K = [2.0, 1.5, 0.5, 0.1]
z = [0.25, 0.25, 0.25, 0.25]
res = rachfordrice(K, z)
β = res.V
y = vapour_compositions(K, z, β)
x = liquid_compositions(K, z, β)
```
"""
function rachfordrice(K, z, β0 = nothing;
                        negative_flash::Bool = false,
                        K_tol= 4*eps(eltype(K)),
                        non_inx = nothing,
                        non_iny = nothing,
                        verbose = nothing)
    # Function to solve Rachdord-Rice mass balance

    RR = RachfordRiceResidual(K,z,non_inx,non_iny)
    RRS = RRState(RR)
    return _rachfordrice(RRS,β0,negative_flash,K_tol,verbose)
end


function rachfordrice(RR::RachfordRiceResidual,β0 = nothing;
                                                negative_flash::Bool = false,
                                                K_tol= 4*eps(eltype(RR)),
                                                verbose = nothing)
    RRS = RRState(RR)
    return _rachfordrice(RRS,β0,negative_flash,K_tol,verbose)
end

function _rachfordrice(RRS::RRState,β0,negative_flash::Bool,K_tol,verbose)
    RR = RRS.RR
    β,status,limits,feasible = rr_initial_state(RRS,β0;K_tol = K_tol)

    _0,_1 = zero(β),one(β)
    nan = _0/_0

    # single-phase shortcuts

    #result is liquid single-phase, (feasible or not), so we return β = 0
    if status == RRLiquid && (!negative_flash)
        return RachfordRiceResult(_0,_1,status,negative_flash,false,feasible,0)

    #result is vapour single-phase, (feasible or not), so we return β = 1
    elseif status == RRVapour && (!negative_flash)
        return RachfordRiceResult(_1,_0,status,negative_flash,false,feasible,0)

    #result is liquid single-phase and no solution is available.
    elseif status == RRLiquid && !feasible
        return RachfordRiceResult(nan,nan,status,negative_flash,false,feasible,0)

    #result is vapour single-phase and no solution is available.
    elseif status == RRVapour && !feasible
        return RachfordRiceResult(nan,nan,status,negative_flash,false,feasible,0)

    #trivial or failure, not feasible otherwise
    elseif status != RREq || !feasible
        return RachfordRiceResult(nan,nan,status,negative_flash,false,feasible,0)
    end

    βv,βl,iter = rrsolve(RR,β,limits) # bracketed Halley when possible

    if !isfinite(βv) || !isfinite(βl)
        if is_verbose(verbose)
            @info "rachford-rice solver: failed to get a result from the iteration" status limits
        end
        return RachfordRiceResult(nan,nan,RRFailure,negative_flash,false,feasible,iter)
    end

    if !negative_flash
        βv = clamp(βv,_0,_1)
        βl = clamp(βl,_0,_1)
    end

    # near-boundary stabilization: if βx is within gate of 0 or 1 and K spans unity,
    # use rr_margin_check to project to boundary.

    if !negative_flash
        β_margin = min(βv, βl)
        gate = cbrt(K_tol)
        if β_margin <= gate
            status_p,βv_margin = rr_margin_check(RR,K_tol,verbose)
            if status_p != RREq
                is_verbose(verbose) && @info "rachfordrice boundary projection" βv βv_margin β_margin gate status_p
                return RachfordRiceResult(βv_margin, 1 - βv_margin,status,negative_flash,true,feasible,iter)
            end
        end
    end

    return RachfordRiceResult(βv,βl,status,negative_flash,true,feasible,iter)
end

function rr_margin_check(RR::RachfordRiceResidual,
                        K_tol= 4*eps(eltype(RR)),
                        verbose = nothing)

    Ktype = eltype(RR)
    _0 = zero(Ktype)
    _1 = one(Ktype)
    F0, Fp0 = residual_first_derivative(RR, zero(Ktype))
    F1, Fp1 = residual_first_derivative(RR, one(Ktype))
    bubble_condition = (abs(F0) <= K_tol) && (F1 < 0) # bubble candidate
    dew_condition = (abs(F1) <= K_tol) && (F0 > 0) # dew candidate

    if is_verbose(verbose)
        δβ = sqrt(eps(Ktype))
        Fδβ = residual(RR,δβ)
        Fp0_num = (Fδβ - F0) / δβ
        @info """ checking boundary conditions:
        F(0)             = $(F0)
        F(1)             = $(F1)
        F'(0)            ≈ $(Fp0_num)
        gate             = $(cbrt(K_tol))
        bubble condition = $(bubble_condition)
        dew condition    = $(dew_condition)
        """
    end

    if bubble_condition && !dew_condition
        β = _0
        status = RRLiquid
        if is_verbose(verbose)
            @info "boundary projection applied: β → 0 (bubble), status = RRLiquid"
        end
    elseif !bubble_condition && dew_condition
        β = _1
        status = RRVapour
        if is_verbose(verbose)
            @info "boundary projection applied: β → 1 (dew), status = RRVapour"
        end
    elseif bubble_condition && dew_condition
        # choose the boundary with smaller linearized step δβ = |F|/|F'|
        δβ0 = abs(F0)/(abs(Fp0) + eps(Ktype))
        δβ1 = abs(F1)/(abs(Fp1) + eps(Ktype))
        if δβ0 <= δβ1
            β = _0
            status = RRLiquid
            if is_verbose(verbose)
                @info "boundary projection via linearized distance: β→0 (δβ0=$(δβ0) ≤ δβ1=$(δβ1))"
            end
        else
            β = _1
            status = RRVapour
            if is_verbose(verbose)
                @info "boundary projection via linearized distance: β→1 (δβ1=$(δβ1) < δβ0=$(δβ0))"
            end
        end
    else
        if is_verbose(verbose)
            @info "boundary projection not conclusive."
        end
        status = RREq
        β = _0/_0
    end
    return status,β
end

#refines a rachford-rice result via Halley iterations
function rrsolve(RR::RachfordRiceResidual,β0,limits = rr_limits(RR,true))

    βmin,βmax = limits
    β_exact,solved = try_rrsolve_23(RR)

    _0,_1 = zero(β_exact),oneunit(β_exact)

    if solved 
        if β_exact < 0.999
            return (β_exact, 1 - β_exact, 0)
        else
            β = β_exact
        end
    else
        β = β0*_1
    end
    
    if !isfinite(β)
        return (β, 1 - β,0)
    end

    rrmin,rrmax = residual(RR,βmin),residual(RR,βmax)

    if sign(rrmin) != sign(rrmax) && βmin <= β <= βmax
        return __rrsolve(RR,β,βmin,βmax)
    end

    rr0 = residual(RR,β)

    if sign(rrmin) != sign(rr0)
        βx = 0.5*(βmin+β)
        βlo,βhi = minmax(β,βmin)
        return __rrsolve(RR,βx,βlo,βhi)
    elseif sign(rrmax) != sign(rr0)
        βx = 0.5*(βmax+β)
        βlo,βhi = minmax(β,βmax)
        return __rrsolve(RR,βx,βlo,βhi)
    else
        return __rrsolve(RR,β)
    end
end

try_rrsolve_23(RR::RachfordRiceResidual) = try_rrsolve_23(RRState(RR))

function try_rrsolve_23(RRS::RRState)
    RR = RRS.RR
    K, z = RR, RR.z
    T = Base.promote_eltype(K, z)

    Z∞ = RRS.Zinf; Z₀ = RRS.Z0
    z1,z2,z3,z4 = T(Z∞),T(Z₀),zero(T),zero(T)
    k1,k2,k3,k4 = T(Inf),zero(T),T(Inf),T(Inf)
    n_max = 2 + iszero(Z₀) + iszero(Z∞)
    n_found = 0
    for i in eachindex(K)
        Ki = K[i]; zi = z[i]
        if zi > 0  && !isinf(Ki) && !iszero(Ki) && !isone(Ki)
            n_found == n_max && return T(NaN),false
            z1,z2,z3,z4 = zi,z1,z2,z3
            k1,k2,k3,k4 = Ki,k1,k2,k3
            n_found += 1
        end
    end
    n_found > n_max && return T(NaN),false
    βmin,βmax = RRS.poles
    degree = n_found + !iszero(Z₀) + !iszero(Z∞)
    a0,a1,a2,a3 = rr_cubic_coeffs((k1,k2,k3,k4),(z1,z2,z3,z4),degree)

    if a3 != 0 && !isnan(a3)
        return cubic_root_in_bracket(a0, a1, a2, a3, βmin, βmax)
    end

    if a2 != 0 && !isnan(a2)
        Δ = a1*a1 - 4*a0*a2
        inva2 = 1/(2*a2)
        Δ2 = sqrt(Δ)*inva2
        x = -a1*inva2
        β1 = x - Δ2
        β2 = x + Δ2
        if βmin < β1 < βmax
            return β1, true
        elseif βmin < β2 < βmax
            return β2, true
        else
            return T(NaN),false
        end
    end

    if a1 != 0 && !isnan(a1)
        β = -a0 / a1
        return (βmin < β < βmax) ? (β, true) : (T(NaN), false)
    end

    return T(NaN),false
end

function rr_cubic_coeffs(K,z,n)
    zz1,zz2,zz3,zz4 = z
    invsumz = 1/sum(z)
    z1,z2,z3,z4 = zz1*invsumz,zz2*invsumz,zz3*invsumz,zz4*invsumz
    T = Base.promote_eltype(K,z)
    k1,k2,k3,k4 = K
    b1,b2,b3,b4 = 1/(1-k1),1/(1-k2),1/(1-k3),1/(1-k4)
    _0 = zero(T)
    if  n == 2
        #0 = a0 + a1β
        a1 = (z1*b2 + z2*b1)
        a0 = z1+z2
        return (T(a0),T(a1),_0,_0)
    elseif n == 3
        #0 = a0 + a1β + a2β^2
        a2 =(z1 + z2 + z3)
        a1 = -b1*(z2 + z3) - b2*(z1 + z3) - b3*(z1 + z2)
        a0 = b1*b2*z3 + b1*b3*z2 + b2*b3*z1
        return (T(a0),T(a1),T(a2),_0)
    elseif n == 4
        #0 = a0 + a1β + a2β^2 + a3β^3
        a3 =(z1 + z2 + z3 + z4)
        a2 = -b1*(z2 + z3 + z4) - b2*(z1 + z3 + z4) - b3*(z1 + z2 + z4) - b4*(z1+z2+z3)
        a1 = b1*b2*(z3+z4) + b1*b3*(z2+z4) + b2*b3*(z1+z4) + b1*b4*(z2+z3) + b2*b4*(z1+z3) + b3*b4*(z1+z2)
        a0 = -b1*b2*b3*z4 - b1*b2*b4*z3 - b1*b3*b4*z2 - b2*b3*b4*z1
        return (T(a0),T(a1),T(a2),T(a3))
    end

    ax = T(NaN)
    return (ax,ax,ax,ax)
end

__rrsolve(RR, βv) =  __rrsolve(RR, βv, oftype(βv,-Inf), oftype(βv,Inf))

function __rrsolve(RR, βv, βmin, βmax)
    T = Base.promote_eltype(RR, βv)
    _0, _1 = zero(T), one(T)


    βl = _1 - βv

    ϵ = eps(T)
    tol = sqrt(ϵ)
    C = isqrt(length(RR)) + 1
    is_vapour = βv <= βl

    if is_vapour
        p, pmin, pmax = βv, βmin, βmax
    else
        p, pmin, pmax = βl, _1 - βmax, _1 - βmin
    end

    iter = 0
    max_iters = 50
    converged = false
    for _ in 1:max_iters
        isfinite(p) || break
        iter += 1
        FO, dFO, d2FO, S = residual_first_second_derivative_error(RR, βv, βl)
        if abs(FO) <= C*ϵ*S
            converged = true
            break
        end

        if (FO < _0) == is_vapour
            pmax = p
        else
            pmin = p
        end

        # Halley step in the primary variable
        den  = 2*dFO*dFO - FO*d2FO
        d = ifelse(is_vapour,1,-1) #isgn
        dp   = -2*FO*dFO*d / den
        pnew = p + dp

        is_halley = pmin < pnew < pmax
        if is_halley
            p = pnew
        elseif isfinite(pmin) && isfinite(pmax)
            pm = (pmin + pmax) / 2
            dp = pm - p
            p = pm
        else
            # bracket is open on the root's side: bisection is meaningless, extrapolate
            dp = -2*FO*d/dFO
            p += dp
        end
        q = _1 - p

        #switch: p > 0.5
        if p > 0.5
            p,q = q,p
            pmin, pmax = _1 - pmax, _1 - pmin
            is_vapour = !is_vapour
        end

        βv, βl = is_vapour ? (p, q) : (q, p)

        if (is_halley && abs(dp) <= tol*abs(p))
            converged = true
            break
        end

        if (pmax - pmin) <= 2*ϵ*abs(p)
            converged = true
            break
        end
    end

    if iter == max_iters || !converged
        βv *= NaN
        βl *= NaN
    end

    return βv, βl, iter
end