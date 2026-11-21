"""
    RachfordRiceStatus

Enumeration describing the phase state found by [`rachfordrice`](@ref).

| value       | meaning |
|:------------|:--------|
| `RRLiquid`  | single-phase liquid: the residual at `β = 0` is negative, so no vapour is formed. |
| `RRVapour`  | single-phase vapour: the residual at `β = 1` is positive, so no liquid is formed. |
| `RREq`      | two-phase (vapour-liquid equilibrium) solution; a root of the Rachford-Rice equation exists. |
| `RRTrivial` | trivial solution: all K-values are equal to one within `K_tol` (`maximum(abs.(log.(K))) <= K_tol`), so the phases are indistinguishable. |
| `RRFailure` | the problem is ill-posed (`NaN` or negative K-values, `NaN` residual) or the iteration did not converge. |

"""
@enum RachfordRiceStatus RRLiquid RRVapour RREq RRFailure RRTrivial

"""
    RachfordRiceResidual(K, z; non_inx = nothing, non_iny = nothing)

Container describing a Rachford-Rice problem: the K-values `K`, the overall molar
composition `z` and the (optional) sets of components that are absent from one of the phases.
The vapour fraction `β` is the root of

    Σᵢ zᵢ (Kᵢ - 1) / (1 + β (Kᵢ - 1)) = 0

(the sum is normalised by `Σᵢ zᵢ`, so `z` does not need to add up to one).

# Arguments
- `K`: vector of K-values `Kᵢ = yᵢ/xᵢ`.
- `z`: vector of overall molar composition, with `length(z) >= length(K)`.

# Keyword arguments
- `non_inx`: `nothing`, or a `Bool` vector marking the components that are not present in the liquid phase. These components are treated as having `Kᵢ = Inf`.
- `non_iny`: `nothing`, or a `Bool` vector marking the components that are not present in the vapour phase. These components are treated as having `Kᵢ = 0`.

A `RachfordRiceResidual` behaves as a read-only iterable container of *effective* K-values.

# Example
```julia
K = [2.0, 1.5, 0.5, 0.1]
z = [0.25, 0.25, 0.25, 0.25]
RR = RachfordRiceResidual(K, z)
res = rachfordrice(RR)
RR[1] #2.0
```
"""
struct RachfordRiceResidual{KK,Z,NX,NY}
    K::KK
    z::Z
    non_inx::NX
    non_iny::NY

    function RachfordRiceResidual(K::KK,z::Z,non_inx::NX,non_iny::NY) where {KK,Z,NX,NY}
        n = length(K)
        checkbounds(z,n)

        if non_inx !== nothing
            checkbounds(non_inx,n)
        end

        if non_iny !== nothing
            checkbounds(non_iny,n)
        end

        return new{KK,Z,NX,NY}(K,z,non_inx,non_iny)
    end
end


const RRSimple{KK,Z} = RachfordRiceResidual{KK,Z,Nothing,Nothing} where {KK,Z}

Base.eltype(RR::RachfordRiceResidual) = eltype(RR.K)
Base.length(RR::RachfordRiceResidual) = length(RR.K)
Base.size(RR::RachfordRiceResidual) = (length(RR.K),)


RachfordRiceResidual(K,z; non_inx = nothing,non_iny = nothing) = RachfordRiceResidual(K,z,non_inx,non_iny)

@inline function Base.getindex(RR::RachfordRiceResidual,i::Int)
    TT = eltype(RR.K)
    return ifelse(rr_any(RR.non_inx,i),TT(Inf),ifelse(rr_any(RR.non_iny,i),TT(0),RR.K[i]))
end

@inline Base.getindex(RR::RRSimple,i::Int) = RR.K[i]
@inline Base.setindex!(RR::RachfordRiceResidual,v,i) = RR.K[i] = v

Base.IteratorSize(::Type{<:RachfordRiceResidual}) = Base.HasShape{1}()

function Base.iterate(RR::RachfordRiceResidual)
    nc = length(RR.K)
    nc == 0 && return nothing
    return iterate(RR,(1,nc))
end

function Base.iterate(RR::RachfordRiceResidual, state)
    i,nc = state
    if i > nc
        return nothing
    else
        return (RR[i],(i + 1,nc))
    end
end

Base.firstindex(RR::RachfordRiceResidual) = firstindex(RR.K)
Base.lastindex(RR::RachfordRiceResidual)  = lastindex(RR.K)
Base.eachindex(RR::RachfordRiceResidual) = eachindex(RR.K)


struct RRState{T,Res<:RachfordRiceResidual}
    RR::Res
    Zinf::T
    Z0::T
    poles::NTuple{2,T}
    sumz::T
    feasible::Int8
end

function RRState(RR::RachfordRiceResidual)
    T = eltype(RR)
    z = RR.z
    sumz = zero(T)
    Kmin,Kmax = one(T)*Inf,-one(T)
    Z0,Zinf = zero(T),zero(T)
    neg_z = false
    for i in eachindex(RR)
        zi,Ki = z[i],RR[i]
        sumz += zi
        if zi < 0
            neg_z = true
        end
        if isinf(Ki)
            Zinf += zi
        elseif iszero(Ki)
            Z0 += zi
        end
        Kmin = min(Ki,Kmin)
        Kmax = max(Ki,Kmax)
    end

    if isnan(Kmin) || isnan(Kmax) || Kmin < 0 || !isfinite(sumz) || neg_z
        feasible = Int8(3)
    elseif Kmin > 1
        feasible = Int8(2)
    elseif Kmax < 1
        feasible = Int8(1)
    else
        feasible = Int8(0)
    end
    _1 = one(T)
    poles = (_1/(_1 - Kmax),_1/(_1 - Kmin))
    return RRState(RR,Zinf,Z0,poles,sumz,feasible)
end

rrstate_absolute_liquid(RRS::RRState) = RRS.feasible == 1
rrstate_absolute_vapour(RRS::RRState) = RRS.feasible == 2
rr_state_kz_domain_err(RRS::RRState) = RRS.feasible == 3

rr_feasible(RRS::RRState) = RRS.feasible == 0
rr_feasible(RR::RachfordRiceResidual) = rr_feasible(RRState(RR))


"""

    rr_feasible(RR::RachfordRiceResidual)::Bool

Returns a `Bool` indicating if the system could have a solution. systems without solution are:

- Systems where all K-values are over 1 (one phase, vapour state)
- Systems where all K-values are below 1 (one phase, liquid state)
- Systems where one or more values of K or z are negative
- Systems where one or more values of K or z are non-finite. K can be `Inf`, but it cannot be `NaN`
"""
rr_feasible

__log_K_norm(RR::RachfordRiceResidual) = __log_K_norm(RR.K,RR.non_inx,RR.non_iny)
__log_K_norm(RR::RRSimple) =maximum(abs ∘ log,RR.K)


"""
    residual(RR::RachfordRiceResidual, βv, βl = 1 - βv)

Evaluate the Rachford-Rice residual function

    F(β) = Σᵢ zᵢ (Kᵢ - 1) / (1 + β (Kᵢ - 1)) / Σᵢ zᵢ

at `β = βv`. The molar vapour fraction of the solution is the root of `F`.

Here `RR` is a [`RachfordRiceResidual`](@ref), `βv` is the molar vapour fraction `β` and `βl`
is the molar liquid fraction, which defaults to `1 - βv`. 

Returns the residual `F(β)`.
"""
function residual(RR::RachfordRiceResidual,βv::B,βl::B = 1 - βv) where B
    K,z = RR,RR.z
    T = Base.promote_eltype(K,z,βv)
    _0,_1 = zero(T),oneunit(T)
    βl_inv = - _1/ βl
    βv_inv = _1/ βv
    res = _0
    res_KD2 = _0
    for i in eachindex(K)
        Ki = K[i]    
        # modification for non-in-y components Ki -> 0
        if iszero(Ki)
            KD2 = -_1
            KD1 = βl_inv - KD2
        # modification for non-in-x components Ki -> ∞
        elseif isinf(Ki)
            KD1 = βv_inv
            KD2 = _0
        # normal Ki
        else
            detKi = inv_fma(Ki,βv,βl)
            KD1 = Ki*detKi
            KD2 = -detKi
        end
        zi = z[i]
        res += zi*KD1
        res_KD2 += zi*KD2
    end
    invsumz = 1/sum(z)
    res += res_KD2
    res *= invsumz
    return res
end

"""
    residual_first_derivative(RR::RachfordRiceResidual, βv, βl = 1 - βv) -> (F, ∂F)

Evaluate the Rachford-Rice residual and its first derivative with respect to `β`:

    F(β) = Σᵢ zᵢ (Kᵢ - 1) / (1 + β (Kᵢ - 1)) / Σᵢ zᵢ
    ∂F(β) = -Σᵢ zᵢ (Kᵢ - 1)² / (1 + β (Kᵢ - 1))² / Σᵢ zᵢ

at `β = βv`. The molar vapour fraction of the solution is the root of `F`.

Here `RR` is a [`RachfordRiceResidual`](@ref), `βv` is the molar vapour fraction `β` and `βl`
is the molar liquid fraction, which defaults to `1 - βv`. 

Returns the tuple `(F, ∂F)`.
"""
function residual_first_derivative(RR::RachfordRiceResidual,βv::B,βl::B = 1 - βv) where B
    K,z = RR,RR.z
    T = Base.promote_eltype(K,z,βv)
    _0,_1 = zero(T),oneunit(T)
    βl_inv = - _1/ βl
    βv_inv = _1/ βv
    res,∂res = _0,_0
    res_KD2 = _0
    for i in eachindex(K)
        Ki = K[i]    
        # modification for non-in-y components Ki -> 0
        if iszero(Ki)
            KD = βl_inv
            KD2 = -_1
            KD1 = βl_inv - KD2
        # modification for non-in-x components Ki -> ∞
        elseif isinf(Ki)
            KD = βv_inv
            KD1 = βv_inv
            KD2 = _0
        # normal Ki
        else
            detKi = inv_fma(Ki,βv,βl)
            KD1 = Ki*detKi
            KD2 = -detKi
            KD = KD1 + KD2
        end
        zi = z[i]
        res += zi*KD1
        KD² = KD*KD
        res_KD2 += zi*KD2
        ∂res -= zi*KD²
    end
    invsumz = 1/sum(z)
    res += res_KD2
    res *= invsumz
    ∂res *= invsumz
    return res,∂res
end

"""
    residual_first_second_derivative(RR::RachfordRiceResidual, βv, βl = 1 - βv) -> (F, ∂F, ∂²F)

Evaluate the Rachford-Rice residual and its first and second derivatives with respect to `β`:

    F(β) = Σᵢ zᵢ (Kᵢ - 1) / (1 + β (Kᵢ - 1)) / Σᵢ zᵢ
    ∂F(β) = -Σᵢ zᵢ (Kᵢ - 1)² / (1 + β (Kᵢ - 1))² / Σᵢ zᵢ
    ∂²F(β) = 2 Σᵢ zᵢ (Kᵢ - 1)³ / (1 + β (Kᵢ - 1))³ / Σᵢ zᵢ

at `β = βv`. The molar vapour fraction of the solution is the root of `F`.

Here `RR` is a [`RachfordRiceResidual`](@ref), `βv` is the molar vapour fraction `β` and `βl`
is the molar liquid fraction, which defaults to `1 - βv`. 

Returns the tuple `(F, ∂F, ∂²F)`
"""
function residual_first_second_derivative(RR::RachfordRiceResidual,βv::B,βl::B = 1 - βv) where B
    res,∂res,∂2res,_ = residual_first_second_derivative_error(RR,βv,βl)
    return res,∂res,∂2res
end

"""
    residual_first_second_derivative_error(RR::RachfordRiceResidual, βv, βl = 1 - βv) -> (F, ∂F, ∂²F, S)

Same as [`residual_first_second_derivative`](@ref), but additionally returns:

    S(β) = Σᵢ zᵢ |Kᵢ - 1| / |1 + β (Kᵢ - 1)| / Σᵢ zᵢ

the sum of the absolute values of the terms that make up `F`. Since `F` is a sum of terms of
both signs, `S` is the natural scale of its floating-point rounding error: a residual with
`|F| <= C * eps * S` cannot be distinguished from zero. The solver uses this as its convergence
criterion.

Here `RR` is a [`RachfordRiceResidual`](@ref), `βv` is the molar vapour fraction `β` and `βl`
is the molar liquid fraction, which defaults to `1 - βv`. 

Returns the tuple `(F, ∂F, ∂²F, S)`.
"""
function residual_first_second_derivative_error(RR::RachfordRiceResidual,βv::B,βl::B = 1 - βv) where B
    K,z = RR,RR.z
    T = Base.promote_eltype(K,z,βv)
    _0,_1 = zero(T),oneunit(T)
    βl_inv = - _1/ βl
    βv_inv = _1/ βv
    res,∂res,∂2res = _0,_0,_0
    res_KD2 = _0
    S = _0
    for i in eachindex(K)
        Ki = K[i] 
        # modification for non-in-y components Ki -> 0
        if iszero(Ki)
            KD = βl_inv
            KD2 = -_1
            KD1 = βl_inv - KD2
        # modification for non-in-x components Ki -> ∞
        elseif isinf(Ki)
            KD = βv_inv
            KD1 = βv_inv
            KD2 = _0
        # normal Ki
        else
            detKi = inv_fma(Ki,βv,βl)
            KD1 = Ki*detKi
            KD2 = -detKi
            KD = KD1 + KD2
        end
        zi = z[i]
        res += zi*KD1
        res_KD2 += zi*KD2
        KD² = KD*KD
        KD³ = KD²*KD
        ∂res -= zi*KD²
        ∂2res += 2*zi*KD³
        S += zi*abs(KD)
    end
    invsumz = 1/sum(z)
    res += res_KD2
    res *= invsumz
    ∂res *= invsumz
    ∂2res *= invsumz
    S *= invsumz
    return res,∂res,∂2res,S
end

function rr_poles(RR::RachfordRiceResidual)
    Kmin,Kmax = extrema(RR)
    return (1/(1 - Kmax), 1/(1 - Kmin))
end

rr_poles(RRS::RRState) = RRS.poles

"""
    βmin,βmax = rr_limits(RR::RachfordRiceResidual, negative_flash::Bool = false)

Whitston-michelsen bounds for the rachford-rice problem.

returns NaN,NaN if the provided K and z values cannot be separated into two phases.
"""
function rr_limits(RR::RachfordRiceResidual, negative_flash::Bool) 
    βmin, βmax = rr_limits(RR)
    if !negative_flash
        _0 = zero(βmin); _1 = one(βmin)
        βmin = max(βmin, _0)
        βmax = min(βmax, _1)
    end
    return βmin, βmax
end

function rr_limits(RRS::RRState, negative_flash::Bool)
    βmin, βmax = rr_limits(RRS)
    if !negative_flash
        _0 = zero(βmin); _1 = one(βmin)
        βmin = max(βmin, _0)
        βmax = min(βmax, _1)
    end
    return βmin, βmax
end

rr_limits(RR::RachfordRiceResidual) = rr_limits(RRState(RR))

function rr_limits(RRS::RRState)
    RR = RRS.RR
    K, z = RR, RR.z
    T = Base.promote_eltype(K, z)
    invsumz = inv(RRS.sumz)

    if !rr_feasible(RRS)
        return T(NaN),T(NaN)
    end
    # βmin_t, βmax_t : tight (composition) bounds from x_i ≤ 1, y_i ≤ 1
    βmin_t = T(-Inf); βmax_t = T(Inf)
    βmin_p,βmax_p = RRS.poles
    for i in eachindex(K)
        Ki = K[i]
        zi = z[i] * invsumz
        if isinf(Ki) || iszero(Ki)                        # pole at 1
            continue
        elseif Ki > 1
            d = Ki - 1
            βmin_t = max(βmin_t, (Ki*zi - 1) / d)
        elseif Ki < 1
            d = 1 - Ki
            βmax_t = min(βmax_t, (1 - zi) / d)
        end
    end
    Z∞,Z0 = RRS.Zinf,RRS.Z0
    # summed tight bounds for K → 0 and K → ∞ dominate the per-component ones
    if Z∞ > 0
        βmin_t = max(βmin_t, Z∞*invsumz)
    end
    if Z0 > 0
        βmax_t = min(βmax_t, 1 - Z0*invsumz)
    end

    # use tight bounds when they form a valid (non-degenerate) bracket;
    # fall back to pole bounds otherwise (e.g. at a single-phase boundary
    # where tight bounds can numerically coincide).
    if βmin_t < βmax_t
        return βmin_t, βmax_t
    else
        return βmin_p, βmax_p
    end
end

"""
    rr_status(RR::RachfordRiceResidual;K_tol = 4*eps(eltype(RR)))

given an [`RachfordRiceResidual`](@ref) return a [RachfordRiceStatus](@ref) enum describing the state of the system.
`K_tol` is used as a measure to see how far apart are the K-values from the trivial solution.
"""
function rr_status(RR::RachfordRiceResidual;K_tol = 4*eps(eltype(RR)))

    g0 = residual(RR,0)
    g1 = residual(RR,1)

    _1 = one(g1)
    _0 = zero(g1)

    if isnan(g0) || isnan(g1)
        status = RRFailure
    elseif __log_K_norm(RR) <= K_tol
        status = RRTrivial
    elseif g0 < 0 #liquid
        status = RRLiquid
    elseif g1 > 0 #gas
        status = RRVapour
    else
        status = RREq
    end
    return status
end

function rr_status(RRS::RRState;K_tol = 4*eps(eltype(RR)))

    if rrstate_absolute_liquid(RRS)
        return RRLiquid
    elseif rrstate_absolute_vapour(RRS)
        return RRVapour
    elseif rr_state_kz_domain_err(RRS)
        return RRFailure
    end
    return rr_status(RRS.RR;K_tol)
end

function rr_initial_state(RRS::RRState,β0 = nothing; K_tol = 4*eps(eltype(RR)))
    
    RR = RRS.RR
    status = rr_status(RRS; K_tol)
    K,z = RR,RR.z
    _1 = one(Base.promote_eltype(K,z))
    _0 = zero(_1)
    nan = _0/_0

    feasible = rr_feasible(RRS)

    if status == RRTrivial || status == RRFailure
        return nan,status,(nan,nan),feasible
    end

    #cannot solve, cannot get bounds
    if !feasible && (status == RRLiquid || status == RRVapour)
        return nan,status,(nan,nan),feasible
    end

    βmin,βmax = rr_limits(RRS)

    #initial point generation
    if β0 !== nothing
        if βmin <= β0 <= βmax && isfinite(β0)
            β = _1*β0
        else
            β = (βmax + βmin)/2
        end
    else
        β = (βmax + βmin)/2
    end

    return β,status,(βmin,βmax),feasible
end
