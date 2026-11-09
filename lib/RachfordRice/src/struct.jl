@enum RachfordRiceStatus RRLiquid RRVapour RREq RRFailure RRTrivial

struct RachfordRiceResidual{KK,Z,NX,NY,T}
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

RachfordRiceResidual(K,z) = RachfordRiceResidual(K,z,nothing,nothing)

@inline function Base.getindex(RR::RachfordRiceResidual,i::Int)
    TT = eltype(RR.K)
    return ifelse(rr_any(RR.non_inx,i),TT(Inf),ifelse(rr_any(RR.non_iny,i),TT(0),RR.K[i]))
end

@inline Base.getindex(RR::RRSimple,i::Int) = RR.K[i]

Base.iterate(RR::RachfordRiceResidual) = iterate(RR, firstindex(RR.K))
Base.iterate(RR::RachfordRiceResidual, state) = iterate(RR, state, RR.K)
function Base.iterate(RR::RachfordRiceResidual, state, K)
    next = iterate(K, state)
    next === nothing && return nothing
    i, newstate = next
    return (RR[i], newstate)
end

Base.firstindex(RR::RachfordRiceResidual) = firstindex(RR.K)
Base.lastindex(RR::RachfordRiceResidual)  = lastindex(RR.K)
Base.eachindex(RR::RachfordRiceResidual) = eachindex(RR.K)

__log_K_norm(RR::RachfordRiceResidual) = __log_K_norm(RR.K,RR.non_inx,RR.non_iny)
__log_K_norm(RR::RRSimple) =maximum(abs ∘ log,RR.K)

function residual(RR::RachfordRiceResidual,β::B) where B
    K,z = RR,RR.z
    T = Base.promote_eltype(K,z,β)
    _0,_1 = zero(T),oneunit(T)
    _0βy = - _1/ (_1 - β)
    _0βx = _1/ β
    res = _0
    res_KD2 = _0
    for i in eachindex(K)
        Ki = K[i]
        #we separate (K-1)/(1-βK) into K/(1-βK) and 1/(1-βK) for numerical reasons.
        #see;
        #Ks_eps_0 = [1.2566703532018493e-21, 3.3506275205339295, 1.0300675710905643e-23, 1.706258568414198e-39, 1.6382855298440747e-20]
        #zs_eps_0 = [0.13754371891028325, 0.29845155687154623, 0.2546683930289046, 0.08177453852283137, 0.22756179266643456]
        if Ki < eps(one(Ki))
            detKi0 = 1 - β
            detKi0 += β*Ki
            detKi = 1/detKi0
        else
            Kim1 = Ki - 1
            detKi0 = 1 + β*Kim1
            detKi = 1/detKi0
        end
        KD1 = Ki*detKi
        KD2 = -detKi
        # modification for non-in-y components Ki -> 0
        if iszero(Ki)
            KD = _0βy
            KD2 = -_1
            KD1 = _0βy - KD2
        end
        # modification for non-in-x components Ki -> ∞
        if isinf(Ki)
            KD = _0βx
            KD1 = _0βx
            KD2 = _0
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

function residual_first_derivative(RR::RachfordRiceResidual,β::B) where B
    K,z = RR,RR.z
    T = Base.promote_eltype(K,z,β)
    _0 = zero(T)
    _0βy = - _1/ (_1 - β)
    _0βx = _1/ β
    res,∂res = _0,_0
    res_KD2 = _0
    for i in eachindex(K)
        Ki = K[i]
        #we separate (K-1)/(1-βK) into K/(1-βK) and 1/(1-βK) for numerical reasons.
        #see;
        #Ks_eps_0 = [1.2566703532018493e-21, 3.3506275205339295, 1.0300675710905643e-23, 1.706258568414198e-39, 1.6382855298440747e-20]
        #zs_eps_0 = [0.13754371891028325, 0.29845155687154623, 0.2546683930289046, 0.08177453852283137, 0.22756179266643456]
        if Ki < eps(one(Ki))
            detKi0 = 1 - β
            detKi0 += β*Ki
            detKi = 1/detKi0
        else
            Kim1 = Ki - 1
            detKi0 = 1 + β*Kim1
            detKi = 1/detKi0
        end
        KD1 = Ki*detKi
        KD2 = -detKi
        KD = KD1 + KD2
        # modification for non-in-y components Ki -> 0
        if iszero(Ki)
            KD = _0βy
            KD2 = -_1
            KD1 = _0βy - KD2
        end
        # modification for non-in-x components Ki -> ∞
        if isinf(Ki)
            KD = _0βx
            KD1 = _0βx
            KD2 = _0
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

function residual_first_second_derivative(RR::RachfordRiceResidual,β::B) where B
    K,z = RR,RR.z
    T = Base.promote_eltype(K,z,β)
    _0 = zero(T)
    _0βy = - _1/ (_1 - β)
    _0βx = _1/ β
    res,∂res,∂2res = _0,_0,_0
    res_KD2 = _0
    for i in eachindex(K)
        Ki = K[i]
        #we separate (K-1)/(1-βK) into K/(1-βK) and 1/(1-βK) for numerical reasons.
        #see;
        #Ks_eps_0 = [1.2566703532018493e-21, 3.3506275205339295, 1.0300675710905643e-23, 1.706258568414198e-39, 1.6382855298440747e-20]
        #zs_eps_0 = [0.13754371891028325, 0.29845155687154623, 0.2546683930289046, 0.08177453852283137, 0.22756179266643456]
        if Ki < eps(one(Ki))
            detKi0 = 1 - β
            detKi0 += β*Ki
            detKi = 1/detKi0
        else
            Kim1 = Ki - 1
            detKi0 = 1 + β*Kim1
            detKi = 1/detKi0
        end
        KD1 = Ki*detKi
        KD2 = -detKi
        KD = KD1 + KD2
        # modification for non-in-y components Ki -> 0
        if iszero(Ki)
            KD = _0βy
            KD2 = -_1
            KD1 = _0βy - KD2
        end
        # modification for non-in-x components Ki -> ∞
        if isinf(Ki)
            KD = _0βx
            KD1 = _0βx
            KD2 = _0
        end
        zi = z[i]
        res += zi*KD1
        res_KD2 += zi*KD2
        KD² = KD*KD
        KD³ = KD²*KD
        ∂res -= zi*KD²
        ∂2res += 2*zi*KD³
    end
    invsumz = 1/sum(z)
    res += res_KD2
    res *= invsumz
    ∂res *= invsumz
    ∂2res *= invsumz
    return res,∂res,∂2res
end

function rr_limits(RR::RachfordRiceResidual)
    K, z = RR, RR.z
    T = Base.promote_eltype(K, z)
    invsumz = inv(sum(z))

    # βmin_t, βmax_t : tight (composition) bounds from x_i ≤ 1, y_i ≤ 1
    # βmin_p, βmax_p : pole bounds, the physical domain of R(β)
    βmin_t = T(-Inf); βmax_t = T(Inf)
    βmin_p = T(-Inf); βmax_p = T(Inf)
    Z∞ = zero(T); Z0 = zero(T)

    for i in eachindex(K)
        Ki = K[i]
        zi = z[i] * invsumz
        if isinf(Ki)
            Z∞ += zi
            βmin_p = max(βmin_p, zero(T))                          # pole at 0
        elseif iszero(Ki)
            Z0 += zi
            βmax_p = min(βmax_p, one(T))                           # pole at 1
        elseif Ki > 1
            d = Ki - 1
            βmin_p = max(βmin_p, -inv(d))                          # 1/(1-Ki)
            βmin_t = max(βmin_t, (Ki*zi - 1) * inv(d))
        elseif Ki < 1
            d = 1 - Ki
            βmax_p = min(βmax_p, inv(d))                           # 1/(1-Ki)
            βmax_t = min(βmax_t, (1 - zi) * inv(d))
        end
    end

    # summed tight bounds for K → 0 and K → ∞ dominate the per-component ones
    if Z∞ > 0
        βmin_t = max(βmin_t, Z∞)
    end
    if Z0 > 0
        βmax_t = min(βmax_t, 1 - Z0)
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

function rr_limits(RR::RachfordRiceResidual, negative_flash::Bool)
    βmin, βmax = rr_limits(RR)
    if !negative_flash
        _0 = zero(βmin); _1 = one(βmin)
        βmin = max(βmin, _0)
        βmax = min(βmax, _1)
    end
    return βmin, βmax
end

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

function rr_initial_state(RR::RachfordRiceResidual,β0 = nothing; K_tol = 4*eps(eltype(RR)))
    status = rr_status(RR; K_tol)
    K,z = RR,RR.z
    _1 = one(Base.promote_eltype(K,z))
    _0 = zero(_1)
    nan = _0/_0

    if status == RRTrivial || status == RRFailure
        return nan,status,(nan,nan)
    elseif status == RRLiquid
        βmin,βmax = rr_limits(RR)
        return _0,status,(βmin,βmax)
    elseif status == RRVapour
        βmin,βmax = rr_limits(RR)
        return _1,status,(βmin,βmax)
    end

    βmin,βmax = rr_limits(RR)

    if β0 !== nothing
        if βmin <= β0 <= βmax
            β = _1*β0
        else
            β = (βmax + βmin)/2
        end
    else
        β = (βmax + βmin)/2
    end

    return β,status,(βmin,βmax)
end
