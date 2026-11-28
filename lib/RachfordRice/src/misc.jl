rr_any(::Nothing) = false
rr_any(rr::AbstractVector) = any(rr)
rr_any(::Nothing, i) = false
rr_any(rr::AbstractVector, i) = rr[i]::Bool

#verbosity: default is nothing, so no side-effects if verbose is not provided.

is_verbose(::Nothing) = false
is_verbose(verbose::Bool) = verbose

function __log_K_norm(K::AbstractVector{T}, non_inx, non_iny) where T
    rr_any(non_inx) && return T(Inf)
    rr_any(non_iny) && return T(Inf)
    return maximum(abs ∘ log, K)
end

"""
    cubic_root_in_bracket(a0, a1, a2, a3, βmin, βmax) -> (β, solved)

Real root of a3·β³ + a2·β² + a1·β + a0 = 0 inside (βmin, βmax).
Assumes a3 ≠ 0. For the RR case a3 ≈ 1 (sum of normalized z's).
"""
function cubic_root_in_bracket(a0, a1, a2, a3, βmin, βmax)
    T = typeof(a0)
    inv_a3 = one(T) / a3

    # Depress: β³ + b·β² + c·β + d = 0, β = t - b/3
    b = a2 * inv_a3
    c = a1 * inv_a3
    d = a0 * inv_a3
    shift = b / 3
    p = c - b*b / 3
    q = 2*b*b*b / 27 - b*c / 3 + d

    # Discriminant of depressed cubic t³ + pt + q = 0
    # D > 0 ⇒ one real root ;  D < 0 ⇒ three real roots
    D = (q/2)^2 + (p/3)^3

    if D >= zero(T)
        # One real root, stable Cardano
        if iszero(q) && iszero(p)
            t = zero(T)
        else
            sqD = sqrt(D)
            w = -q/2 + copysign(sqD, -q/2)   # larger of the two
            u = cbrt(w)
            v = -p / (3*u)
            t = u + v
        end
        β = t - shift
        return (βmin < β < βmax) ? (β, true) : (zero(T)/zero(T), false)
    else
        # Three real roots, trigonometric
        r = sqrt(-p/3)
        cos_arg = clamp(-q/(2*r*r*r), -one(T), one(T))
        θ = acos(cos_arg) / 3
        two_pi_over_3 = 2 * (π * one(T)) / 3
        for k in 0:2
            t = 2 * r * cos(θ - k * two_pi_over_3)
            β = t - shift
            βmin < β < βmax && return (β, true)
        end
        return (zero(T)/zero(T), false)
    end
end

function inv_fma(x, y, z)
    r = fma(x, y, z)
    return 1/r
end

macro public(sym)
    sym isa Symbol || error("@public expects a single symbol, got: $sym")
    if VERSION >= v"1.11.0-DEV.469"
        return esc(Expr(:public, sym))
    else
        return nothing
    end
end

function discriminant(a::T, b::T, c::T) where {T<:AbstractFloat}
    p  = b * b
    q  = a * c
    dp = fma(b, b, -p)   # exact rounding error of p
    dq = fma(a, c, -q)   # exact rounding error of q
    return (p - q) + (dp - dq)
end

discriminant(a, b, c) = b*b - a*c
