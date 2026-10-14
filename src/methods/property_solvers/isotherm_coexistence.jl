#=

------------- Coexistence functions --------------------------

Their main use is to bracket a range where there exists a maxwell loop.

=#

const BracketState{T} = NamedTuple{(:v,:p,:dpdv,:d2pdv2),NTuple{4,T}}

@inline function d2pdv2_state(model,xx,T,z,::Type{TT}) where TT
    _f, _df, _d2f = p∂p∂2p(model, xx, T, z)
    return BracketState{TT}((TT(xx),TT(_f),TT(_df),TT(_d2f)))
end

@noinline function coexistence_volumes_bracketed_bracket_log(msg::String,data)
    @info msg data
    return nothing
end

function coexistence_volumes_bracketed(model, T, z; kwargs...)
    lbv = lb_volume(model,T,z)
    vl = 1.001*lbv
    vhi = 1e5*lbv
    return coexistence_volumes_bracketed(model,T,z,vl,vhi;kwargs...)
end

function coexistence_volumes_bracketed(model, T, z, _v_lb, _v_ub, vm0 = nothing;
        tol_v    = 1e-12,
        max_iters::Int  = 100,
        branches::Symbol = :all,
        strikes::Int = 3,
        verbose::Bool  = false)

    TT = promote_type(Base.promote_eltype(model,T,z),Base.promote_eltype(_v_lb,_v_ub))
    v_lb,v_ub = minmax(TT(_v_lb),TT(_v_ub))
    #function barrier: the implementation is specialized on TT.
    return _coexistence_volumes_bracketed(model,T,z,v_lb,v_ub,vm0,TT,
                                          tol_v,max_iters,branches,strikes,verbose)
end

function _coexistence_volumes_bracketed(model,T,z,v_lb::TT,v_ub::TT,vm0,::Type{TT},
        tol_v,max_iters::Int,branches::Symbol,strikes::Int,verbose::Bool)::Tuple{BracketState{TT},BracketState{TT},BracketState{TT},Symbol} where TT

    ∑z = TT(sum(z))
    nan = TT(NaN)
    s0 = BracketState{TT}((nan,nan,nan,nan))
    sl,sv = d2pdv2_state(model,v_lb,T,z,TT),d2pdv2_state(model,v_ub,T,z,TT)

    _bracket_hermite(sa,sb) = Solvers.hermite5_poly(sa.v, sb.v, sa.p, sb.p, sa.dpdv, sb.dpdv, sa.d2pdv2, sb.d2pdv2)
    d2pdrho2(s,∑z) = s.v * s.v * s.v *(s.v * s.d2pdv2 + 2 * s.dpdv) / (∑z*∑z)

    verbose && coexistence_volumes_bracketed_bracket_log(lazy"Clapeyron.coexistence_volumes_bracketed: start",
        (T=T, v_lb=sl.v, p_lb=sl.p, dpdv_lb=sl.dpdv, v_ub=sv.v, p_ub=sv.p, dpdv_ub=sv.dpdv, vm0=vm0, branches=branches))
    #non finite pressures (e.g. volume below the co-volume), or non-positive vapour pressure.
    if !isfinite(sl.p) || !isfinite(sv.p) || sv.p <= 0 || sl.v < 0 || sv.v < 0
        return sl,sv,s0,:failure
    end

    # sl is unstable, sv is stable
    if sl.dpdv > 0 && sv.dpdv < 0
        verbose && coexistence_volumes_bracketed_bracket_log(lazy"Clapeyron.coexistence_volumes_bracketed: case sl unstable, sv stable",(;))
        is_vapour(branches) && (return s0,sv,sl,:vapour)
        #NOTE (type stability): the TT(...) conversion keeps vl_new concrete.
        vl_new = TT(coexistence_expand_volume_bracket(model,T,z,sl.v,:liquid))
        verbose && coexistence_volumes_bracketed_bracket_log(lazy"Clapeyron.coexistence_volumes_bracketed: expanded liquid bracket",(v_old=sl.v, v_new=vl_new))
        isnan(vl_new) && (return sl,sv,sl,:failure)
        sl_new = d2pdv2_state(model,vl_new,T,z,TT)
        return sl_new,sv,sl,:eq
    end

    # sl is stable, sv is unstable
    if sv.dpdv > 0 && sl.dpdv < 0
        verbose && coexistence_volumes_bracketed_bracket_log(lazy"Clapeyron.coexistence_volumes_bracketed: case sl stable, sv unstable",(;))
        is_liquid(branches) && (return sl,s0,sv,:liquid)
        vv_new = TT(coexistence_expand_volume_bracket(model,T,z,sv.v,:vapour))
        verbose && coexistence_volumes_bracketed_bracket_log(lazy"Clapeyron.coexistence_volumes_bracketed: expanded vapour bracket",(v_old=sv.v, v_new=vv_new))
        isnan(vv_new) && (return sl,sv,sv,:failure)
        sv_new = d2pdv2_state(model,vv_new,T,z,TT)
        return sl,sv_new,sv,:eq
    end

    # sl is unstable, sv is unstable (??)
    if sv.dpdv > 0 && sl.dpdv > 0
        verbose && coexistence_volumes_bracketed_bracket_log(lazy"Clapeyron.coexistence_volumes_bracketed: case sl and sv unstable",(;))
        vv_new = TT(coexistence_expand_volume_bracket(model,T,z,sv.v,:vapour))
        vl_new = TT(coexistence_expand_volume_bracket(model,T,z,sl.v,:liquid))
        verbose && coexistence_volumes_bracketed_bracket_log(lazy"Clapeyron.coexistence_volumes_bracketed: expanded both brackets",
            (vl_old=sl.v, vl_new=vl_new, vv_old=sv.v, vv_new=vv_new))
        (isnan(vv_new) || isnan(vl_new)) && (return sl,sv,s0,:failure)
        sv_new = d2pdv2_state(model,vv_new,T,z,TT)
        sl_new = d2pdv2_state(model,vl_new,T,z,TT)
        #pick the unstable point closer to the inflection of dpdv (d2pdv2 = 0).
        #v*d2pdv2/dpdv is dimensionless, so both sides are comparable.
        smid = ifelse(abs(sl.v*sl.d2pdv2/sl.dpdv) < abs(sv.v*sv.d2pdv2/sv.dpdv),sl,sv)
        return sl_new,sv_new,smid,:eq
    end

    if vm0 === nothing
        sm = s0
    else
        #NOTE (type stability): vm0 is converted to TT before use.
        vm0_TT = TT(vm0)
        if !isnan(vm0_TT)
            sm = d2pdv2_state(model,vm0_TT,T,z,TT)
            sm.dpdv > 0 && (return sl,sv,sm,:eq)
            d2pdrho2_m0 = d2pdrho2(sm,∑z)
            if d2pdrho2_m0 > zero(TT)
                sl = sm
            else
                sv = sm
            end
        else
            sm = s0
        end
    end
    strike_count = 0
    status = :failure
    niter = 0
    vapour_updated = false #true if the last update was on the vapour side. if so, the next middle point is generated in density variables.
    verbose && coexistence_volumes_bracketed_bracket_log(lazy"Clapeyron.coexistence_volumes_bracketed: entering hermite loop",
        (v_l=sl.v, v_v=sv.v, dpdv_l=sl.dpdv, dpdv_v=sv.dpdv))

    for i in 1:max_iters
        niter = i

        lnv_converged,_ = Solvers.convergence(log(sl.v),log(sv.v),zero(TT),tol_v)
        lnv_converged && break
        d2pdρ2l = d2pdrho2(sl,∑z)
        d2pdρ2v = d2pdrho2(sv,∑z)
        #NOTE (type stability): vm is declared as TT, so it is concrete in both branches.
        local vm::TT
        if vapour_updated
            #density formulation: the vapour is the lower density bound, the liquid the upper one.
            #we interpolate q = -p, so that the unstable region has dqdρ > 0, with a maximum of dqdρ
            #at the inflexion point, the same sign convention that _find_vm expects for volumes.
            ρl,ρv = ∑z/sl.v,∑z/sv.v
            dqdρl,dqdρv = sl.dpdv*sl.v*sl.v/∑z,sv.dpdv*sv.v*sv.v/∑z
            d2qdρ2l = -d2pdρ2l
            d2qdρ2v = -d2pdρ2v
            poly = Solvers.hermite5_poly(ρv, ρl, -sv.p, -sl.p, dqdρv, dqdρl, d2qdρ2v, d2qdρ2l)
            dpoly = Solvers.polyder(poly)
            ρm = TT(_find_vm(dpoly, ρv, ρl, -ρv)) #polynomial is in (ρ - ρv)
            vm = ∑z/ρm #NaN if ρm is NaN, vm < 0 if ρm < 0
        else
            # Hermite interpolant on the current bounds, volume formulation
            poly = _bracket_hermite(sl,sv)
            c = - sl.v #translation of the hermite poly. v(eos) + c = v(hermite)
            dpoly = Solvers.polyder(poly)

            #get vm = d2pdv2 = 0 via hermite polynomials
            vm = TT(_find_vm(dpoly, sl.v, sv.v,c))
        end

        #of hermite interpolant failed, fall back to bisection.
        #if _find_vm failed, then the interpolant failed to generate a proper maxwell loop
        #if this happens too much times, then, at this temperature, maybe there is no maxwell loop at all
        #for that, we use a strike system. 4 strikes and you are out
        hermite_failed = isnan(vm) || vm < 0
        if hermite_failed
            strike_count += 1
            vm = sqrt(sl.v*sv.v)
        end

        verbose && coexistence_volumes_bracketed_bracket_log(lazy"Clapeyron.coexistence_volumes_bracketed: iteration",
            (i=i, v_l=sl.v, v_v=sv.v, vm=vm, hermite_failed=hermite_failed, strikes=strikes))

        strike_count >= strikes && break

        sm = d2pdv2_state(model,vm,T,z,TT)

        #return if we found our middle point.
        if sm.dpdv > 0
            status = :eq
            verbose && coexistence_volumes_bracketed_bracket_log(lazy"Clapeyron.coexistence_volumes_bracketed: found unstable middle point",
                (i=i, vm=sm.v, dpdv=sm.dpdv, p=sm.p))
            break
        end

        #update bounds
        #Anderson and Prausnitz criteria to distinguish between liquid and gas, good enough here
        #for an stricter criteria, pip was used, but it requires evaluating the EoS again.
        d2pdrho2_m = d2pdrho2(sm,∑z)

        verbose && coexistence_volumes_bracketed_bracket_log(lazy"Clapeyron.coexistence_volumes_bracketed: stable middle point",
            (i=i, vm=sm.v, dpdv=sm.dpdv, d2pdrho2=d2pdrho2_m, moves=(d2pdrho2_m > zero(TT) ? :sl : :sv)))

        if d2pdrho2_m > zero(TT)
            sl = sm
        else
            if !isfinite(sm.p) || sm.p <= zero(TT)
                status = :failure
                break
            end
            sv = sm
            vapour_updated = true
            continue
        end
        vapour_updated = false
    end

    verbose && coexistence_volumes_bracketed_bracket_log(lazy"Clapeyron.coexistence_volumes_bracketed: end",
        (status=status, iterations=niter, strikes=strikes, v_l=sl.v, v_v=sv.v, v_m=sm.v))

    return sl,sv,sm,status
end

function _find_vm(dpoly, lo::K, hi::K,c::K = zero(K)) where K
   #dpoly : dpdv
    d2poly = Solvers.polyder(dpoly) #: d2pdv2
    d3poly = Solvers.polyder(d2poly)#: d3pdv3

    lb = lo + c
    ub = hi + c

    nr, v1, v2, v3 = Solvers.real_roots3(d2poly)

    if (lb <= v1 <= ub) && evalpoly(v1, dpoly) > 0 && evalpoly(v1, d3poly) < 0
        return v1 - c
    elseif nr > 1 && (lb <= v3 <= ub) && evalpoly(v3, dpoly) > 0 && evalpoly(v3, d3poly) < 0
        return v3 - c
    elseif nr > 2 && (lb <= v2 <= ub) && evalpoly(v2, dpoly) > 0 && evalpoly(v2, d3poly) < 0
        return v2 - c
    else
        return zero(K) / zero(K)
    end
end

function coexistence_expand_volume_bracket(model,T,z,v,phase)

    #generate a bracket point that actually has dpdv < 0

    TT = Base.promote_eltype(model,T,z,v)
    vx = TT(v)
    rhox = 1/vx
    is_gas = is_vapour(phase)
    npts = 10
    if is_gas
        xrange = range(rhox,zero(rhox),npts + 1)
    else
        lbv = TT(lb_volume(model,T,z))
        xrange = range(log2(vx),log2(lbv),npts + 1)
    end

    for i in 2:npts
        if is_gas
            vxi = 1/xrange[i]
        else
            vxi = exp2(xrange[i])
        end
        _,dpdv_i = p∂p∂V(model,vxi,T,z)
        dpdv_i < 0 && return vxi
    end
    return zero(vx)/zero(vx)
end

#=
Returns:
bracket_l = (sl_outer,sl) :two states that bracket the coexistence in the dense ("liquid") branch of the isotherm
bracket_v = (sv,sv_outer) : two states that bracket the coexistence in the dilute ("vapor") branch of the isotherm
middle   =  (sml,smv)     : unstable states, but (sl,sml) brackets the liquid spinodal and (smv,sv) brackets the vapour spinodal
spinodal_p = (pl,pv)      : spinodal pressures (estimated from the hermite interpolation, or true spinodals if the fallback was used).

if an isogibbs equilibria is done on the specified tempereature T and composition z (Clapeyron.edge_pressure), then bracket_l brackets the liquid pressures and volumes
and bracket_v brackets the vapour compositions and volumes.

On failure (no maxwell loop, middle point not unstable, failed spinodal or outer bracket), all states are returned as given/partially refined and the spinodal pressures are NaN.

The function is split in small pieces, each one is a top level function without closures (better inference and compilation times):
 - _good_bracket_virial_seed: tighter initial vapour state from the second virial coefficient.
 - _good_bracket_inner: step 1, refinement of the inner brackets (hermite cycle + newton fallback).
   - _good_bracket_step: one hermite iteration on one side: takes a stable and an unstable state, generates a new point,
     replaces the bound of the same stability and checks if the side is done. the same function is used for the liquid (volume
     coordinates to generate the point) and the vapour (density coordinates to generate the point).
   - _good_bracket_close: end of one side: newton fallback if degenerate, final spinodal pressure otherwise.
 - _good_bracket_outer: step 2, outer brackets (liquid extension, zero pressure approximation, vapour ideal gas point and post-check).
=#

function coexistence_get_good_bracket(model,T,z;max_iters_hermite = 15,rtol_spinodal = 1e-8)
    s0l,s0v,s0m,status =  coexistence_volumes_bracketed(model,T,z)
    return coexistence_get_good_bracket(model,T,z,s0l,s0v,s0m,max_iters_hermite,rtol_spinodal)
end

function coexistence_get_good_bracket(model,T,z,s0l,s0v,s0m,max_iters_hermite = 15,rtol_spinodal = 1e-8)
    TT = Base.promote_eltype(model,T,z,s0l.v,s0v.v,s0m.v)
    ∑z = TT(sum(z))
    RT = Rgas(model)*T

    #stable inner brackets
    sl = (v = TT(s0l.v),p = TT(s0l.p),dpdv = TT(s0l.dpdv))
    sv = (v = TT(s0v.v),p = TT(s0v.p),dpdv = TT(s0v.dpdv))

    #unstable bounds
    sm = (v = TT(s0m.v),p = TT(s0m.p),dpdv = TT(s0m.dpdv))
    sml,smv = sm,sm


    #the middle point must be an unstable state, otherwise there is nothing to refine.
    if !(sm.dpdv > 0)
        return (return (sl,sl),(sv,sv),(sm,sm),(TT(NaN),TT(NaN)))
    end

    #stable outer brackets: they are saved before the virial seed, those must keep the original sv.
    sl_outer,sv_outer = sl,sv

    #Step 0: virial seed for the vapour side
    #the second virial coefficient is estimated from the stable vapour input state: p*v = nRT*(1 + B/v_molar).
    #the truncated virial EoS has its spinodal at v_molar = -2B, a good guess for the vapour spinodal.
    #the point is evaluated and placed by the sign of dpdv, so a bad estimate of B costs at most one evaluation.
    B_est = (sv.v/∑z)*(sv.p*sv.v/(∑z*RT) - 1)
    v_vir = -2*B_est*∑z
    if isfinite(B_est) && B_est < 0 && sm.v < v_vir < sv.v
        svir = dpdv_state(model,v_vir,T,z)
        if isfinite(svir.p) && svir.dpdv < 0
            sv = svir #virial point is now the best vapour inner bound
        elseif isfinite(svir.p) && svir.dpdv > 0 && svir.v > sm.v
            smv = svir #virial point is now the best vapour unstable bound
        end
        #virial test failed.
    end

    #Step 1: inner brackets, sl.p < sv.p, (sl,sml) brackets the liquid spinodal and (smv,sv) the vapour spinodal
    sl,sml,smv,sv,p̄l,p̄v,ok = _good_bracket_inner(model,T,z,sl,sml,smv,sv,max_iters_hermite,rtol_spinodal)
    ok || (return (sl,sl),(sv,sv),(sm,sm),(TT(NaN),TT(NaN)))

    #Step 2: outer brackets, such as sl_outer.p > p̄v, sv_outer.p < p̄l (*)
    # on the case of a liquid bracket that reaches negative pressures, we proceed to use the zero pressure approximation.
    sl_outer,sv_outer,ok = _good_bracket_outer(model,T,z,sl,sl_outer,sv_outer,p̄l,p̄v)
    ok || (return (sl,sl),(sv,sv),(sm,sm),(TT(NaN),TT(NaN)))

    return (sl_outer,sl),(sv,sv_outer),(sml,smv),(p̄l,p̄v)
end

function dpdv_state(model,v,T,z)
    TT = Base.promote_eltype(model,T,z,v)
    p,dpdv = p∂p∂V(model,v,T,z)
    return (v = TT(v),p = TT(p),dpdv = TT(dpdv))
end

function _hermite_spinodal(xa,xb,pa,pb,dpa,dpb)
    #xa < xb, in any coordinate x (v or ρ). one state is stable and the other unstable.
    #returns (xs,ps): the extremum of the hermite-3 interpolant of p(x) between them
    #(an estimate of the spinodal), in the original coordinate. (NaN,NaN) if there is no root inside.
    poly = Solvers.hermite3_poly(xa,xb,pa,pb,dpa,dpb) #polynomial in (x - xa)
    dpoly = Solvers.polyder(poly)
    _,x1,x2 = Solvers.real_roots2(dpoly)
    hi = xb - xa
    xs = ifelse(zero(hi) <= x1 <= hi,x1,x2)
    if zero(hi) <= xs <= hi
        return xs + xa,evalpoly(xs,poly)
    end
    return zero(hi)/zero(hi),zero(hi)/zero(hi)
end

function _good_bracket_step(model,T,z,sst,sun,p̄,∑z,rtol_spinodal,density::Bool)
    #one iteration on one side (liquid or vapour). takes a stable and an unstable state, generates a new point near the spinodal,
    #replaces the bound of the same stability and checks if we are done with this side.
    #returns the new stable and unstable states, the spinodal pressure estimate, the newton seed, and the done/fallback flags.

    if density
        #work in density formulation
        dpdρ_sst,dpdρ_sun = -sst.dpdv*sst.v*sst.v/∑z,-sun.dpdv*sun.v*sun.v/∑z
        ρs,ps = _hermite_spinodal(∑z/sst.v,∑z/sun.v,sst.p,sun.p,dpdρ_sst,dpdρ_sun) #ρst < ρun
        vs = ∑z/ρs
    else
        #work in volume formulation
        vs,ps = _hermite_spinodal(sst.v,sun.v,sst.p,sun.p,sst.dpdv,sun.dpdv) #vst < vun
    end

    if !isnan(vs)
        p̄ = ps
    else #no valid root: bisect in log(v) (same as in log(ρ)) and keep the previous estimate of p̄
        vs = sqrt(sst.v*sun.v)
    end

    snew = dpdv_state(model,vs,T,z)
    new_is_unstable = snew.dpdv > 0
    vprev = new_is_unstable ? sun.v : sst.v #bound that is about to be replaced

    #convergence: log(v) == -log(ρ), sqrt(ρ1*ρ2) == 1/sqrt(v1*v2), so it does not depend on the coordinate.
    #if the new point is unstable, the stable bound does not move, so we compare it with the middle of the new bracket.
    vx1 = sst.v #old stable bound
    if new_is_unstable
        vx2 = sqrt(vx1*snew.v)
        sun = snew
    else
        vx2 = snew.v
        sst = snew
    end
    done,fallback = false,false
    converged,_ = Solvers.convergence(log(vx1),log(vx2),zero(vx2),1e-12)
    if converged
        done = true
    elseif abs(snew.v - vprev) <= 4*eps(snew.v)
        #stagnation: the new point reproduces the bound it replaces. if the bracket is already tight,
        #we reached the pressure resolution and we accept it. otherwise, the cycle is degenerate.
        done = true
        fallback = abs(sun.v - sst.v) > rtol_spinodal*max(sun.v,sst.v)
    end
    return sst,sun,p̄,vs,done,fallback
end

function _good_bracket_close(model,T,z,sst,sun,p̄,vseed,fallback::Bool,∑z,density::Bool,rtol_spinodal)
    #end of one side. degenerate sides use the true spinodal (safeguarded newton on dpdv inside the bracket).
    #the other ones recompute the spinodal pressure from the final bracket: by construction of the hermite interpolant,
    #p̄l <= sl.p and sv.p <= p̄v. returns the final bracket, the spinodal pressure and a success flag.
    if fallback
        sst,sun,sn = pure_spinodal_newton_bracket(model,T,z,sst,sun,vseed;rtol = rtol_spinodal)
        isnan(sn.v) && return sst,sun,p̄,false
        p̄ = sn.p
    else
        if density
            #work in density formulation
            dpdρ_sst,dpdρ_sun = -sst.dpdv*sst.v*sst.v/∑z,-sun.dpdv*sun.v*sun.v/∑z
            _,p̄f = _hermite_spinodal(∑z/sst.v,∑z/sun.v,sst.p,sun.p,dpdρ_sst,dpdρ_sun) #ρst < ρun
        else
            #work in volume formulation
            _,p̄f = _hermite_spinodal(sst.v,sun.v,sst.p,sun.p,sst.dpdv,sun.dpdv) #vst < vun
        end

        p̄ = isnan(p̄f) ? p̄ : p̄f
    end
    return sst,sun,p̄,true
end

function _good_bracket_inner(model,T,z,sl,sml,smv,sv,max_iters_hermite,rtol_spinodal)
    ∑z = sum(z)
    TT = Base.promote_eltype(model,T,z,sl.v,sml.v,smv.v,∑z)
    #estimations of the spinodal pressures
    p̄l,p̄v = TT(NaN),TT(NaN)

    done_l = false #sl.p < 0 is checked after the first iteration, so p̄l is always estimated.
    done_v = false
    cross = false #we stop once we find crossing volumes and the inner bounds cross the middle spinodal pressure.

    #the hermite cycle is cheap but can degenerate (one side stuck, or no convergence in max_iters_hermite).
    #in that case, that side falls back to the true spinodal (safeguarded newton on dpdv, more expensive).
    fallback_l,fallback_v = false,false
    vseed_l,vseed_v = TT(NaN),TT(NaN) #last hermite estimate of the spinodal volumes, used as newton seeds.

    #find crossing inner bounds, such as sl.p < sv.p. sl also stops if sl.p < 0
    for i in 1:max_iters_hermite
        #do at least one iteration to get the spinodal estimation
        cross && done_l && done_v && i > 1 && break

        if !done_l
            sl,sml,p̄l,vseed_l,dl,fl = _good_bracket_step(model,T,z,sl,sml,p̄l,∑z,rtol_spinodal,false)
            done_l |= dl
            fallback_l |= fl
        end

        if sl.p < 0
            done_l = true
        end

        if !done_v
            sv,smv,p̄v,vseed_v,dv,fv = _good_bracket_step(model,T,z,sv,smv,p̄v,∑z,rtol_spinodal,true)
            done_v |= dv
            fallback_v |= fv
        end

        if sl.p < sv.p
            cross = true
        end

        if cross
            #liquid side: sl.p < 0 (above), unclipped middle pressure (only if p̄l >= 0), or convergence (in the step).
            #vapour side: only convergence of (smv,sv).
            #if p̄l < 0 the midpoint rule is skipped: step 2 needs sl.p < 0 to use the zero pressure approximation,
            #and that is already covered by the sl.p < 0 check above.
            p̄m = (p̄l + p̄v)/2
            p̄l >= 0 && sl.p < p̄m && (done_l = true)
        end
    end

    #sides that did not converge in the cheap cycle are degenerate.
    fallback_l |= !done_l
    fallback_v |= !done_v


    sl,sml,p̄l,ok_l = _good_bracket_close(model,T,z,sl,sml,p̄l,vseed_l,fallback_l,∑z,false,rtol_spinodal)
    ok_l || return sl,sml,smv,sv,p̄l,p̄v,false
    sv,smv,p̄v,ok_v = _good_bracket_close(model,T,z,sv,smv,p̄v,vseed_v,fallback_v,∑z,true,rtol_spinodal)
    ok_v || return sl,sml,smv,sv,p̄l,p̄v,false

    #the spinodal pressures bound the stable states: p̄l <= sl.p and sv.p <= p̄v. this protects against a failed root or roundoff.
    p̄l = isnan(p̄l) ? sl.p : min(p̄l,sl.p)
    p̄v = isnan(p̄v) ? sv.p : max(p̄v,sv.p)


    #the inner bounds must cross (sl.p < sv.p), otherwise there is no maxwell loop at this temperature.
    ok = sl.p < sv.p
    return sl,sml,smv,sv,p̄l,p̄v,ok
end

function _good_bracket_outer(model,T,z,sl,sl_outer,sv_outer,p̄l,p̄v)
    #step 2: outer brackets, such as sl_outer.p > p̄v >= psat >= p̄l_mod > sv_outer.p
    #returns the outer liquid and vapour states and a success flag.
    nRT = sum(z)*Rgas(model)*T

    #--- liquid: push sl_outer above p̄v ---
    if sl_outer.p < p̄v
        for i in 1:100
            #p̄v is the pressure of a converged spinodal, so it is already an upper bound of psat. no margin needed:
            sl_outer.p > p̄v && break
            vᵢ,pᵢ,dpdvᵢ = sl_outer.v,sl_outer.p,sl_outer.dpdv
            #p(v) is convex in the liquid branch, so this step lands at or above the target.
            Δᵢ = (p̄v-pᵢ)/(vᵢ*dpdvᵢ)
            sl_outer = dpdv_state(model,vᵢ*exp(Δᵢ),T,z) #logVᵢ = logVᵢ + Δᵢ
        end
    end

    #--- vapour: pressure target p̄l_mod, a lower bound of psat ---
    #we trust the liquid bracket more than p̄l: if the bracket reaches negative pressures, we are at low temperatures and
    #the zero pressure approximation is used: the liquid volume is the one at p = 0 and the vapour is an ideal gas (lnϕ = 0).
    #then, lnϕ(liquid) = 0 is solved for p. p̄l_mod is a lower bound of psat (the liquid fugacity grows with p, the vapour one is lower than p if Z < 1)
    if sl.p < 0
        if !(sl_outer.p > 0)
            #"Clapeyron.coexistence_get_good_bracket: outer liquid does not reach positive pressures, returning failure"
            return sl_outer,sv_outer,false
        end
        vl0 = coexistence_interpolate_volume(sl_outer,sl,zero(sl.p))
        ares_l0 = a_res(model, vl0, T, z)
        lnf_l0 = ares_l0 - 1 + log(nRT/vl0)
        p̄l_mod = exp(lnf_l0)
    elseif p̄l > 0
        p̄l_mod = p̄l
    else
        # "Clapeyron.coexistence_get_good_bracket: spinodal pressure is not positive but the liquid bracket does not reach negative pressures, returning failure"
        return sl_outer,sv_outer,false
    end

    if sv_outer.p > p̄l_mod
        #just use ideal gas approximation, after all, with the other bracket we can generate a good EoS in density
        sv_outer = dpdv_state(model,nRT/p̄l_mod,T,z)
    end

    #post-check: the outer vapour point needs p < p̄l_mod (<= psat). with an ideal gas volume this holds if Z < 1 at that state,
    #which is not guaranteed for every EoS. on the stable dilute branch p decreases with v, so we expand the volume until it holds.
    for i in 1:10
        sv_outer.p < p̄l_mod && break
        sv_outer = dpdv_state(model,2*sv_outer.v,T,z)
    end
    if !(sv_outer.p < p̄l_mod)
        return sl_outer,sv_outer,false
    end
    return sl_outer,sv_outer,true
end

function pure_spinodal_newton_bracket(model, T, z, sst, sun, vs0 = nothing; rtol = 1e-8, max_iters = 25)
    #sst is a stable state (dpdv < 0) and sun an unstable one (dpdv > 0), in any order of volume.
    #returns (sst,sun,s): the final bracket (still one stable and one unstable state) and
    #s, the best spinodal state (NaN on failure).
    stable_is_lo = sst.v < sun.v
    sa,sb = stable_is_lo ? (sst,sun) : (sun,sst) #sa.v < sb.v
    nan = zero(sa.v)/zero(sa.v)
    snan = (v = nan, p = nan, dpdv = nan)
    sa.dpdv*sb.dpdv > 0 && return sst,sun,snan
    fscale = max(abs(sa.dpdv),abs(sb.dpdv))
    vs = (vs0 !== nothing && sa.v < vs0 < sb.v) ? vs0 : (sa.v + sb.v)/2
    sbest = snan
    for j in 1:max_iters
        pj,fs,d2pj = p∂p∂2p(model,vs,T,z)
        s = (v = vs, p = pj, dpdv = fs)
        #update the bracket first, so the sign change stays inside
        #(the sign of sa.dpdv never changes, so sa keeps its stability).
        if fs*sa.dpdv <= 0
            sb = s
        else
            sa = s
        end
        Δ = fs/d2pj
        vn = vs - Δ
        if abs(fs) <= rtol*fscale || abs(sb.v - sa.v) <= rtol*sb.v || abs(Δ) <= 4*eps(vs)
            sbest = s
            break
        end
        #safeguarded step: newton only if it stays strictly inside the bracket (also catches NaN/Inf)
        vs = (sa.v < vn < sb.v) ? vn : (sa.v + sb.v)/2
    end

    if stable_is_lo
        return sa,sb,sbest
    else
        return sb,sa,sbest
    end
end

function coexistence_interpolate_volume(sa, sb, p)
    slo, shi = sa.p > sb.p ? (sb, sa) : (sa, sb)
    dp = shi.p - slo.p
    βp = (p - slo.p) / dp
    if !(zero(βp) <= βp <= one(βp))
        throw(DomainError(lazy"Clapeyron.coexistence_interpolate_volume: pressure $(p) outside bounds ($(slo.p),$(shi.p))"))
    end

    # --- endpoint classification ----------------------------------------
    # dpdv < 0 on both stable branches. Near a spinodal dpdv → 0⁻, so the
    # "soft" (spinodal-adjacent) endpoint is the one with smaller |dpdv|
    #dlo, dhi = abs(slo.dpdv), abs(shi.dpdv)
    #soft_is_lo = dlo < dhi

    # --- Hermite data in x = ln v ---------------------------------------
    x0, x1 = log(slo.v), log(shi.v)
    m0 = inv(slo.v * slo.dpdv)   # dx/dp at slo
    m1 = inv(shi.v * shi.dpdv)   # dx/dp at shi
    Δ  = (x1 - x0) / dp          # secant slope in x

    # --- Fritsch-Carlson monotonicity test ------------------------------
    # α, β ≥ 0 and α² + β² ≤ 9 ⇒ cubic in x(p) is provably monotone on [slo.p, shi.p].
    # The sign test matters if one endpoint is not stable (dpdv > 0): α² + β² ≤ 9 alone does not
    # imply monotonicity in that case. If it holds, direct evaluation is safe and cheap.
    α = m0 / Δ
    β = m1 / Δ
    fc_ok = α >= 0 && β >= 0 && (α*α + β*β) <= 9

    # --- fast path: direct cubic in x(p) --------------------------------
    if fc_ok
        poly = Solvers.hermite3_poly(slo.p, shi.p, x0, x1, m0, m1)
        return exp(evalpoly(p - slo.p, poly))
    end

    # --- slow path: invert p(v) via cubic root solve --------------------
    # p(v) is monotone on this interval (both dpdv < 0), so the cubic has
    # at most one root in [shi.v, slo.v] and root_in_bounds is reliable.
    poly = Solvers.hermite3_poly(shi.v, slo.v, shi.p, slo.p, shi.dpdv, slo.dpdv)
    k0, k1, k2, k3 = poly
    poly2 = (k0 - p, k1, k2, k3)
    res = Solvers.real_roots3(poly2)
    cv = -shi.v
    in_bounds, vx = Solvers.root_in_bounds(res, shi.v + cv, slo.v + cv)
    if in_bounds
        return vx - cv
    end

    # --- final fallback: linear in ln v ---------------------------------
    # Only reachable if the cubic root finder missed, which should be rare.
    return exp(x0 + βp*(x1 - x0))
end