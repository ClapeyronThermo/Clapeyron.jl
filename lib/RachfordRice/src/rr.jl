function rachfordrice(K, z, β0 = nothing;
                        negative_flash::Bool = false
                        K_tol= 4*eps(eltype(RR)), 
                        non_inx = nothing, 
                        non_iny = nothing, 
                        verbose = false)
    # Function to solve Rachdord-Rice mass balance

    RR = RachfordRiceResidual(K,z,non_inx,non_iny)

    return _rachfordrice(RR,β0,negative_flash,Ktol,verbose)
end


function rachfordrice(RR::RachfordRiceResidual,β0 = nothing;
                                                negative_flash::Bool = false
                                                K_tol= 4*eps(eltype(RR)),
                                                verbose = false)

    return _rachfordrice(RR,β0,negative_flash,Ktol,verbose)
end

function _rachfordrice(RR::RachfordRiceResidual,β0,negative_flash::Bool,K_tol,verbose)
    
    β,status,limits = rr_initial_state(RR,β0;K_tol = K_tol)
    _0,_1 = zero(β),one(β)
    # single-phase shortcuts
    
    if status == RRLiquid
        return _0
    elseif status == RRVapour
        return _1
    elseif status != RREq
        return _0/_0
    end

    
    βx = rrsolve(RR,β,limits,negative_flash) # bracketed Halley when possible
    
    if !negative_flash
        βx = clamp(βx,_0,_1)
    end

    # near-boundary stabilization: if βx is within gate of 0 or 1 and K spans unity,
    # use rr_margin_check to project to boundary.
    Kmin,Kmax = extrema(RR)
    _11 = one(Kmin)
    if (Kmin < _11 && Kmax > _11)
        β_margin = min(βx, _1 - βx)
        gate = cbrt(K_tol)
        if β_margin <= gate
            status_p,βb = rr_margin_check(RR;K_tol,verbose)
            if status_p != RREq
                verbose && @info "rachfordrice boundary projection" βx βb β_margin gate status_p
                return clamp(βb,_0,_1)
            end
        end
    end
    return βx
end

function rr_margin_check(RR::RachfordRiceResidual,
                        K_tol= 4*eps(eltype(RR)), 
                        verbose = false)

    Ktype = eltype(RR)
    _0 = zero(Ktype)
    _1 = one(Ktype)
    F0 = residual(RR,0) # F(0)
    F1 = residual(RR,1) # F(1)
    cond0 = (abs(F0) <= K_tol) && (F1 < 0) # bubble candidate
    cond1 = (abs(F1) <= K_tol) && (F0 > 0) # dew candidate

    if is_verbose(verbose)
        δβ = sqrt(eps(Ktype))
        Fδβ = residual(RR,δβ)
        Fp0_num = (Fδβ - F0) / δβ
        @info """ checking boundary conditions:
        F(0)             = $(F0)
        F(1)             = $(F1)
        F'(0)            ≈ $(Fp0_num)
        gate             = $(cbrt(K_tol))
        bubble condition = $(cond0)
        dew condition    = $(cond1)
        """
    end

    if cond0 && !cond1
        β = _0
        status = RRLiquid
        if is_verbose(verbose)
            @info "boundary projection applied: β→0 (bubble), status=RRLiquid"
        end
    elseif !cond0 && cond1
        β = _1
        status = RRVapour
        if is_verbose(verbose)
            @info "boundary projection applied: β→1 (dew), status=RRVapour"
        end
    elseif cond0 && cond1
        # choose the boundary with smaller linearized step δβ = |F|/|F'|
        Fp0,Fp1 = _0,_0
        @inbounds for i in eachindex(K,z)
            Δ   = K[i] - _1
            zi  = z[i]
            Fp0 -= zi*(Δ*Δ)
            Fp1 -= zi*(Δ*Δ)/(K[i]*K[i])
        end
        δβ0 = abs(F0)/(abs(Fp0)+eps(Ktype))
        δβ1 = abs(F1)/(abs(Fp1)+eps(Ktype))
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
function rrsolve(RR::RachfordRiceResidual,β0,limits = rr_limits(RR,true),negative_flash::Bool = false)
    
    β_exact,status_exact = try_rrsolve_23(RR)
    if status_exact == RREq && return β_exact
        if negative_flash
            return β_exact
        else
            b0,b1 = zero(β_exact), oneunit(β_exact)
            return clamp(β_exact,b0,b1)
        end
    end

    βmin,βmax = limits

    if !negative_flash
        βmin = max(zero(βmin),βmin)
        βmax = min(zero(βmax),βmax)
    end

    _0 = zero(Base.promote_eltype(RR,β0,βmin,βmax))
    _1 = one(_0)

    
    β = β0
    
    if !isfinite(β)
        return β
    end

    rrmin,rrmax = residual(RR,βmin),residual(RR,βmax)
    F = Base.Fix1(residual_first_second_derivative,RR)

    if rrmin*rrmax <= 0 && βmin <= β <= βmax
        return guarded_halley(F,β,βmin,βmax)
    end

    rr0 = residual(RR,β)

    if rrmin*rr0 <= 0
        βx = 0.5*(βmin+β)
        βlo,βhi = minmax(β,βmin)
        return guarded_halley(F,βx,βlo,βhi)
    elseif rrmax*rr0 <= 0
        βx = 0.5*(βmax+β)
        βlo,βhi = minmax(β,βmax)
        return guarded_halley(F,βx,βlo,βhi)
    else
        return guarded_halley(F,β)
    end
end

@inline function trysolve23_phys(β::T, d1, d2, d3, c, Z∞, Z₀) where T
    Z∞ > zero(T) && β <= zero(T) && return false
    Z₀ > zero(T) && β >= one(T)  && return false
    c >= 1 && one(T) + β*d1 <= zero(T) && return false
    c >= 2 && one(T) + β*d2 <= zero(T) && return false
    c >= 3 && one(T) + β*d3 <= zero(T) && return false
    return true
end

function try_rrsolve_23(RR::RachfordRiceResidual)
    K, z = RR, RR.z
    T = zero(Base.promote_eltype(K, z))
    nc = length(K)

    # ---------- pass 1: classify ----------
    Z∞ = zero(T); Z₀ = zero(T)
    n∞ = 0; n₀ = 0; n₁ = 0
    has_light = false
    has_heavy = false
    for i in eachindex(K)
        Ki = K[i]; zi = z[i]
        (Ki < 0 || isnan(Ki)) && return (zero(T), RRFailure)
        if zi > 0
            if isinf(Ki)
                n∞ += 1; Z∞ += zi
            elseif iszero(Ki)
                n₀ += 1; Z₀ += zi
            elseif isone(Ki)
                n₁ += 1
            elseif Ki > 1
                has_light = true
            else
                has_heavy = true
            end
        else
            n₁ += 1
        end
    end

    degP = nc - 1 - n₁ - max(n₀ - 1, 0) - max(n∞ - 1, 0)
    degP > 2 && return (zero(T), RRFailure)

    L_finite = (Z∞ > zero(T)) || has_light
    U_finite = (Z₀ > zero(T)) || has_heavy

    if !L_finite && !U_finite
        return (zero(T), RRTrivial)
    elseif !U_finite
        return (one(T), RRVapour)
    elseif !L_finite
        return (zero(T), RRLiquid)
    end

    # ---------- pass 2: collect ≤ 3 significant (d, z) pairs ----------
    d1 = zero(T); d2 = zero(T); d3 = zero(T)
    z1 = zero(T); z2 = zero(T); z3 = zero(T)
    c = 0
    for i in eachindex(K)
        Ki = K[i]; zi = z[i]
        if zi > 0 && isfinite(Ki) && !iszero(Ki) && !isone(Ki)
            c += 1
            d = Ki - one(T)
            if     c == 1; d1 = d; z1 = zi
            elseif c == 2; d2 = d; z2 = zi
            else;          d3 = d; z3 = zi
            end
        end
    end

    # ---------- build P(β) = p0 + p1 β + p2 β² ----------
    p0 = zero(T); p1 = zero(T); p2 = zero(T)

    if Z∞ > zero(T) && Z₀ > zero(T)
        if c == 0
            p0 = Z∞; p1 = -(Z∞ + Z₀)
        else # c == 1
            p0 = Z∞
            p1 = d1 * (Z∞ + z1) - (Z∞ + Z₀)
            p2 = -d1 * (Z∞ + Z₀ + z1)
        end
    elseif Z∞ > zero(T)
        if c == 0
            p0 = Z∞
        elseif c == 1
            p0 = Z∞; p1 = d1 * (Z∞ + z1)
        else # c == 2
            S0 = z1*d1 + z2*d2
            p0 = Z∞
            p1 = Z∞*(d1 + d2) + S0
            p2 = d1 * d2 * (Z∞ + z1 + z2)
        end
    elseif Z₀ > zero(T)
        if c == 0
            p0 = -Z₀
        elseif c == 1
            p0 = z1*d1 - Z₀
            p1 = -d1 * (Z₀ + z1)
        else # c == 2
            S0 = z1*d1 + z2*d2
            S1 = d1 * d2 * (z1 + z2)
            p0 = S0 - Z₀
            p1 = (S1 - S0) - Z₀*(d1 + d2)
            p2 = -Z₀*d1*d2 - S1
        end
    else
        if c == 1
            p0 = z1*d1
        elseif c == 2
            p0 = z1*d1 + z2*d2
            p1 = d1 * d2 * (z1 + z2)
        else # c == 3
            p0 = z1*d1 + z2*d2 + z3*d3
            p1 = z1*d1*(d2 + d3) + z2*d2*(d1 + d3) + z3*d3*(d1 + d2)
            p2 = d1 * d2 * d3 * (z1 + z2 + z3)
        end
    end

    # ---------- solve, pick physical root ----------
    if iszero(p2)
        iszero(p1) && return (zero(T), RRFailure)
        β = -p0 / p1
        return (β, trysolve23_phys(β, d1, d2, d3, c, Z∞, Z₀) ? RREq : RRFailure)
    end

    two  = one(T) + one(T)
    disc = p1*p1 - two*two*p2*p0
    disc < zero(T) && return (zero(T), RRFailure)
    sq = sqrt(disc)
    q  = -(p1 + copysign(sq, p1)) / two

    if iszero(q)  # p0 = p1 = 0: double root at β = 0
        return trysolve23_phys(zero(T), d1, d2, d3, c, Z∞, Z₀) ?
               (zero(T), RREq) : (zero(T), RRFailure)
    end

    β1 = q / p2
    β2 = p0 / q
    trysolve23_phys(β1, d1, d2, d3, c, Z∞, Z₀) && return (β1, RREq)
    trysolve23_phys(β2, d1, d2, d3, c, Z∞, Z₀) && return (β2, RREq)
    return (zero(T), RRFailure)
end