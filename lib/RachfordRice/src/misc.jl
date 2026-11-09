rr_any(::Nothing) = false
rr_any(rr::AbstractVector) = any(rr)
rr_any(::Nothing,i) = false
rr_any(rr::AbstractVector,i) = rr[i]::Bool

#verbosity: default is nothing, so no side-effects if verbose is not provided.

is_verbose(::Nothing) = false
is_verbose(verbose::Bool) = verbose

function __log_K_norm(K::AbstractVector{T},non_inx,non_iny) where T
    rr_any(non_inx) && return T(Inf)
    rr_any(non_iny) && return T(Inf)
    return maximum(abs ∘ log,K)
end

function guarded_halley(f, x0::T) where T
    guarded_halley(f, x0, T(-Inf), T(Inf))
end

function guarded_halley(f, x0::T, xlo::T, xhi::T) where T
    zero_T = zero(T); one_T = one(T); two_T = one_T + one_T
    ft = eps(T)
    bracket = isfinite(xlo) && isfinite(xhi) && xlo < xhi

    lo, hi, flo = xlo, xhi, zero_T
    if bracket
        flo = first(f(lo))
        fhi = first(f(hi))
        iszero(flo) && return lo
        iszero(fhi) && return hi
        (flo < zero_T) == (fhi < zero_T) &&
            throw(ArgumentError("f(xlo) and f(xhi) must have opposite signs"))
        xt = eps(T) * two_T * (abs(lo) + abs(hi) + one_T)
        x = clamp(x0, lo, hi)
    else
        xt = zero_T
        x = x0
    end

    for _ in 1:100
        fx, dfx, d2fx = f(x)
        abs(fx) <= ft && return x
        bracket && (hi - lo) <= xt && return (lo + hi) / two_T

        if bracket
            if (fx < zero_T) == (flo < zero_T)
                lo = x; flo = fx
            else
                hi = x
            end
        end

        denom = two_T * dfx * dfx - fx * d2fx
        xn = (iszero(denom) || !isfinite(denom)) ?
             (bracket ? (lo + hi) / two_T : return(x)) :
             x - two_T * fx * dfx / denom

        if !isfinite(xn) || xn == x || (bracket && !(lo < xn < hi))
            bracket || return x
            xn = (lo + hi) / two_T
        end
        x = xn
    end

    return x
end