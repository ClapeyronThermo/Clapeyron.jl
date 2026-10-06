#=

------------- Coexistence functions --------------------------

Their main use is to bracket a range where there exists a maxwell loop.

=#

function coexistence_volumes_bracketed(model, T, z, _v_lb, _v_ub, vm0 = nothing;
        tol_v    = 1e-12,
        maxiter  = 100,
        branches = :all,
        verbose  = false)

    TT = Base.promote_eltype(model,T,z,_v_lb,_v_ub)
    v_lb,v_ub = minmax(TT(_v_lb),TT(_v_ub))
    ∑z = TT(sum(z))

    function p_state(xx)
        _f, _df, _d2f = p∂p∂2p(model, xx, T, z)
        return (v = TT(xx),p = TT(_f),dpdv = TT(_df),d2pdv2 = TT(_d2f))
    end

    hermite_from_state(sa,sb) = Solvers.hermite5_poly(sa.v, sb.v, sa.p, sb.p, sa.dpdv, sb.dpdv, sa.d2pdv2, sb.d2pdv2)

    sl,sv = p_state(v_lb),p_state(v_ub)
    s0 = (v = TT(NaN),p = TT(NaN),dpdv = TT(NaN),d2pdv2 = TT(NaN))

    verbose && @info "coexistence_volumes_bracketed: start" T=T v_lb=sl.v p_lb=sl.p dpdv_lb=sl.dpdv v_ub=sv.v p_ub=sv.p dpdv_ub=sv.dpdv vm0=vm0 branches=branches

    #non finite pressures (e.g. volume below the co-volume), or non-positive vapour pressure.
    if !isfinite(sl.p) || !isfinite(sv.p) || sv.p <= 0
        return sl,sv,s0,:failure
    end

    # sl is unstable, sv is stable
    if sl.dpdv > 0 && sv.dpdv < 0
        verbose && @info "coexistence_volumes_bracketed: case sl unstable, sv stable"
        is_vapour(branches) && (return s0,sv,sl,:vapour)
        vl_new = coexistence_expand_volume_bracket(model,T,z,sl.v,:liquid)
        verbose && @info "coexistence_volumes_bracketed: expanded liquid bracket" v_old=sl.v v_new=vl_new
        isnan(vl_new) && (return sl,sv,sl,:failure)
        sl_new = p_state(vl_new)
        return sl_new,sv,sl,:eq
    end

    # sl is stable, sv is unstable
    if sv.dpdv > 0 && sl.dpdv < 0
        verbose && @info "coexistence_volumes_bracketed: case sl stable, sv unstable"
        is_liquid(branches) && (return sl,s0,sv,:liquid)
        vv_new = coexistence_expand_volume_bracket(model,T,z,sv.v,:vapour)
        verbose && @info "coexistence_volumes_bracketed: expanded vapour bracket" v_old=sv.v v_new=vv_new
        isnan(vv_new) && (return sl,sv,sv,:failure)
        sv_new = p_state(vv_new)
        return sl,sv_new,sv,:eq
    end

    # sl is unstable, sv is unstable (??)
    if sv.dpdv > 0 && sl.dpdv > 0
        verbose && @info "coexistence_volumes_bracketed: case sl and sv unstable"
        vv_new = coexistence_expand_volume_bracket(model,T,z,sv.v,:vapour)
        vl_new = coexistence_expand_volume_bracket(model,T,z,sl.v,:liquid)
        verbose && @info "coexistence_volumes_bracketed: expanded both brackets" vl_old=sl.v vl_new=vl_new vv_old=sv.v vv_new=vv_new
        (isnan(vv_new) || isnan(vl_new)) && (return sl,sv,s0,:failure)
        sv_new = p_state(vv_new)
        sl_new = p_state(vl_new)
        #pick the unstable point closer to the inflection of dpdv (d2pdv2 = 0).
        #v*d2pdv2/dpdv is dimensionless, so both sides are comparable.
        smid = ifelse(abs(sl.v*sl.d2pdv2/sl.dpdv) < abs(sv.v*sv.d2pdv2/sv.dpdv),sl,sv)
        return sl_new,sv_new,smid,:eq
    end

    if vm0 === nothing
        sm = s0
    else
        if !isnan(vm0)
            sm = p_state(TT(vm0))
            sm.dpdv > 0 && (return sl,sv,sm,:eq)
            d2pdrho2_m0 = sm.v * sm.v * sm.v *(sm.v * sm.d2pdv2 + 2 * sm.dpdv) / (∑z*∑z)
            if d2pdrho2_m0 > zero(TT)
                sl = sm
            else
                sv = sm
            end 
        else
            sm = s0
        end
    end

    strikes = 0
    status = :failure
    niter = 0

    verbose && @info "coexistence_volumes_bracketed: entering hermite loop" v_l=sl.v v_v=sv.v dpdv_l=sl.dpdv dpdv_v=sv.dpdv

    for i in 1:maxiter
        niter = i

        lnv_converged,_ = Solvers.convergence(log(sl.v),log(sv.v),zero(TT),tol_v)
        lnv_converged && break

        # Hermite interpolant on the current bounds
        poly = hermite_from_state(sl,sv)
        c = - sl.v #translation of the hermite poly. v(eos) + c = v(hermite)
        dpoly = Solvers.polyder(poly)

        #get vm = d2pdv2 = 0 via hermite polynomials
        vm = _find_vm(dpoly, sl.v, sv.v,c)

        #of hermite interpolant failed, fall back to bisection.
        #if _find_vm failed, then the interpolant failed to generate a proper maxwell loop
        #if this happens too much times, then, at this temperature, maybe there is no maxwell loop at all
        #for that, we use a strike system. 4 strikes and you are out
        hermite_failed = isnan(vm)
        if isnan(vm)
            strikes += 1
            vm = sqrt(sl.v*sv.v)
        end

        verbose && @info "coexistence_volumes_bracketed: iteration" i=i v_l=sl.v v_v=sv.v vm=vm hermite_failed=hermite_failed strikes=strikes

        strikes >= 4 && break

        sm = p_state(vm)

        #return if we found our middle point.
        if sm.dpdv > 0
            status = :eq
            verbose && @info "coexistence_volumes_bracketed: found unstable middle point" i=i vm=sm.v dpdv=sm.dpdv p=sm.p
            break
        end

        #update bounds
        #Anderson and Prausnitz criteria to distinguish between liquid and gas, good enough here
        #for an stricter criteria, pip was used, but it requires evaluating the EoS again.
        d2pdrho2_m = sm.v * sm.v * sm.v *(sm.v * sm.d2pdv2 + 2 * sm.dpdv) / (∑z*∑z)

        verbose && @info "coexistence_volumes_bracketed: stable middle point" i=i vm=sm.v dpdv=sm.dpdv d2pdrho2=d2pdrho2_m moves=(d2pdrho2_m > zero(TT) ? :sl : :sv)

        if d2pdrho2_m > zero(TT)
            sl = sm
        else
            sv = sm
        end
    end

    verbose && @info "coexistence_volumes_bracketed: end" status=status iterations=niter strikes=strikes v_l=sl.v v_v=sv.v v_m=sm.v

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

function coexistence_get_good_bracket(model,T,z,s0l,s0v,s0m;verbose = false)
    #=
    Returns:
    bracket_l = (sl_outer,sl) :two states that bracket the coexistence in the dense ("liquid") branch of the isotherm
    bracket_v = (sv,sv_outer) : two states that bracket the coexistence in the dilute ("vapor") branch of the isotherm
    middle   =  (sml,smv)     : unstable states, but (sl,sml) brackets the liquid spinodal and (smv,sv) brackets the vapour spinodal
    spinodal_p = (pl,pv)      : extrapolated spinodal pressures.

    if an isogibbs equilibria is done on the specified tempereature T and composition z (Clapeyron.edge_pressure), then bracket_l brackets the liquid pressures and volumes
    and bracket_v brackets the vapour compositions and volumes.
    =#
    TT = Base.promote_eltype(model,T,z,s0l.v,s0v.v,s0m.v)
    ∑z = TT(sum(z))

    #density coordinates for the vapour side: ρ = ∑z/v, dp/dρ = -dpdv*v²/∑z
    dpdρ(s) = -s.dpdv*s.v*s.v/∑z

    function p_state(xx)
        _f, _df = p∂p∂V(model, xx, T, z)
        return (v = TT(xx),p = TT(_f),dpdv = TT(_df))
    end

    #hermite_from_state_old(sa,sb) = Solvers.hermite5_poly(sa.v, sb.v, sa.p, sb.p, sa.dpdv, sb.dpdv, sa.d2pdv2, sb.d2pdv2)
    hermite_from_state(sa,sb) = Solvers.hermite3_poly(sa.v, sb.v, sa.p, sb.p, sa.dpdv, sb.dpdv)

    #stable inner brackets
    sl = (v = s0l.v,p = s0l.p,dpdv = s0l.dpdv)
    sv = (v = s0v.v,p = s0v.p,dpdv = s0v.dpdv)

    #metastable bounds
    sm = (v = s0m.v,p = s0m.p,dpdv = s0m.dpdv)
    sml,smv = sm,sm

    verbose && @info "coexistence_get_good_bracket: start" T=T sum_z=∑z sl=sl sv=sv sm=sm

    #the middle point must be an unstable state, otherwise there is nothing to refine.
    if !(sm.dpdv > 0)
        verbose && @info "coexistence_get_good_bracket: middle point is not unstable, returning early" dpdv_m=sm.dpdv
        return (sl,sl),(sv,sv),(sm,sm),(TT(NaN),TT(NaN))
    end

    #stable outer brackets
    sl_outer,sv_outer = sl,sv

    #estimations of the spinodal pressures
    p̄l,p̄v = TT(NaN),TT(NaN)

    done_l = false #sl.p < 0 is checked after the first iteration, so p̄l is always estimated.
    done_v = false
    cross = false #we stop once we find crossing volumes and the inner bounds cross the middle spinodal pressure.
    niter = 0
    #Step 1: find crossing inner bounds, such as sl.p < sv.p. sl also stops if sl.p < 0
    for i in 1:60
        #do at least one iteration to get the spinodal estimation
        cross && done_l && done_v && i > 1 && break
        niter = i
        if !done_l
            poly_l = hermite_from_state(sl,sml) #hermite-3 interpolation to find spindodal
            dpoly_l = Solvers.polyder(poly_l) #dpdv polynomial is a quadratic, the root inside the bounds is the spindodal prediction
            _,vl1,vl2 = Solvers.real_roots2(dpoly_l)
            cl = -sl.v #translation from hermite implementation.
            vsl = ifelse((sl.v + cl) <= vl1 <= (sml.v + cl),vl1,vl2)
            root_ok_l = (sl.v + cl) <= vsl <= (sml.v + cl)
            if root_ok_l
                p̄l = evalpoly(vsl,poly_l)
            else #no valid root: bisect in log(v) and keep the previous estimate of p̄l
                vsl = sqrt(sl.v*sml.v) + cl
            end
            ssl = p_state(vsl - cl)
            verbose && @info "get_good_bracket liquid" i=i v_sl=sl.v v_sml=sml.v p_sl=sl.p dpdv_sl=sl.dpdv dpdv_sml=sml.dpdv root_candidates=(vl1 - cl,vl2 - cl) root_ok=root_ok_l v_spinodal_est=vsl - cl pbar_l=p̄l v_new=ssl.v p_new=ssl.p dpdv_new=ssl.dpdv moves=(ssl.dpdv > 0 ? :sml : :sl)
            if ssl.dpdv > 0
                vxl1 = sl.v
                vxl2 = sqrt(sl.v*ssl.v)
                sml = ssl
            else
                vxl1 = sl.v
                vxl2 = ssl.v
                sl = ssl

            end
            
            lnvl_converged,_ = Solvers.convergence(log(vxl1),log(vxl2),zero(TT),1e-12)

            if lnvl_converged
                done_l = true
            end
        end

        if sl.p < 0
            done_l = true
        end

        if !done_v
            #vapour side works in density: p(ρ) is concave from the ideal gas limit up to the
            #unstable inflexion, so the spinodal is a plain maximum (no change of sign of d2p/dρ2,
            #unlike d2p/dv2 which changes sign between the spinodal and the stable inflexion).
            ρsv,ρsmv = ∑z/sv.v,∑z/smv.v #ρsv < ρsmv
            ρs,ps = _hermite_spinodal(ρsv,ρsmv,sv.p,smv.p,dpdρ(sv),dpdρ(smv))
            root_ok_v = !isnan(ρs)
            if isnan(ρs) #no valid root: bisect in log(ρ) and keep the previous estimate of p̄v
                ρs = sqrt(ρsv*ρsmv)
            else
                p̄v = ps
            end
            ssv = p_state(∑z/ρs)
            ρssv = ∑z/ssv.v
            verbose && @info "get_good_bracket vapour" i=i rho_sv=ρsv rho_smv=ρsmv v_sv=sv.v v_smv=smv.v p_sv=sv.p dpdrho_sv=dpdρ(sv) dpdrho_smv=dpdρ(smv) root_ok=root_ok_v rho_spinodal_est=ρs pbar_v=p̄v v_new=ssv.v p_new=ssv.p dpdv_new=ssv.dpdv moves=(ssv.dpdv > 0 ? :smv : :sv)
            if ssv.dpdv > 0
                ρxv1 = ρsv
                ρxv2 = sqrt(ρsv*ρssv)
                smv = ssv
            else
                ρxv1 = ρsv
                ρxv2 = ρssv
                sv = ssv
            end
            lnρv_converged,_ = Solvers.convergence(log(ρxv1),log(ρxv2),zero(TT),1e-12)

            if lnρv_converged
                done_v = true
            end
        end

        if sl.p < sv.p
            cross = true
        end
        
        if cross
            #liquid side: sl.p < 0 (above), unclipped middle pressure, or convergence (above).
            #vapour side: only convergence of (smv,sv) in log(ρ).
            p̄m = (p̄l + p̄v)/2
            sl.p < p̄m && (done_l = true)
            verbose && @info "get_good_bracket middle pressure" i=i pbar_m=p̄m sl_p=sl.p liquid_below_middle=(sl.p < p̄m) sv_p=sv.p
        end
        verbose && @info "get_good_bracket iteration summary" i=i cross=cross done_l=done_l done_v=done_v sl_v=sl.v sml_v=sml.v smv_v=smv.v sv_v=sv.v sl_p=sl.p sv_p=sv.p pbar_l=p̄l pbar_v=p̄v
    end

    
    verbose && @info "get_good_bracket: step 1 finished" iterations=niter cross=cross done_l=done_l done_v=done_v

    #final spinodal estimates, computed from the final brackets (sl,sml) and (smv,sv).
    #by construction of the hermite interpolant, p̄l <= sl.p and sv.p <= p̄v. the min/max
    #only protects against a failed root (fallback) or roundoff.
    _,p̄l_f = _hermite_spinodal(sl.v,sml.v,sl.p,sml.p,sl.dpdv,sml.dpdv) #volume coordinates
    _,p̄v_f = _hermite_spinodal(∑z/sv.v,∑z/smv.v,sv.p,smv.p,dpdρ(sv),dpdρ(smv)) #density coordinates
    p̄l = isnan(p̄l_f) ? p̄l : p̄l_f
    p̄v = isnan(p̄v_f) ? p̄v : p̄v_f
    p̄l = isnan(p̄l) ? sl.p : min(p̄l,sl.p)
    p̄v = isnan(p̄v) ? sv.p : max(p̄v,sv.p)

    verbose && @info "get_good_bracket: final spinodal pressure estimates" pbar_l=p̄l pbar_v=p̄v sl_p=sl.p sv_p=sv.p

    #Step 2: find crossing outer bounds, such as sl.p > p̄v, sv.p < p̄l (*)
    # on the case of negative p̄l, we proceed to use the zero volume approximation.

    verbose && @info "get_good_bracket: step 2 liquid outer" v_outer=sl_outer.v p_outer=sl_outer.p pbar_v=p̄v needs_extension=(sl_outer.p < p̄v)
    if sl_outer.p < p̄v
        for i in 1:100
            sl_outer.p > p̄v && break
            vᵢ,pᵢ,dpdvᵢ = sl_outer.v,sl_outer.p,sl_outer.dpdv
            Δᵢ = (1.01p̄v-pᵢ)/(vᵢ*dpdvᵢ) #1.01 is to allow passing p̄v without traps.
            sl_outer = p_state(vᵢ*exp(Δᵢ)) #logVᵢ = logVᵢ + Δᵢ
            verbose && @info "get_good_bracket: liquid outer extension" i=i delta_lnv=Δᵢ v_outer=sl_outer.v p_outer=sl_outer.p dpdv_outer=sl_outer.dpdv
        end
    end

    RT = Rgas(model)*T
    nRT = ∑z*RT
    if sl_outer.p < 0 || p̄l < 0
        #zero_pressure approximation.
        vl0 = coexistence_interpolate_volume(sl_outer,sl,zero(TT),true)
        ares_l0 = a_res(model, vl0, T, z)
        lnϕ_l0 = ares_l0 - 1 + log(nRT/vl0)
        p̄l_mod = exp(lnϕ_l0)
        verbose && @info "get_good_bracket: zero pressure approximation" vl0=vl0 ares_l0=ares_l0 pbar_l_mod=p̄l_mod pbar_l=p̄l
    else
        p̄l_mod = p̄l
    end

    verbose && @info "get_good_bracket: step 2 vapour outer" v_outer=sv_outer.v p_outer=sv_outer.p pbar_l_mod=p̄l_mod needs_extension=(sv_outer.p > p̄l_mod)
    if sv_outer.p > p̄l_mod
        #just use ideal gas approximation, after all, with the other bracket we can generate a good EoS in density
        sv_outer = p_state(nRT/p̄l_mod)
        verbose && @info "get_good_bracket: ideal gas vapour outer" v_outer=sv_outer.v p_outer=sv_outer.p dpdv_outer=sv_outer.dpdv
    end

    verbose && @info "get_good_bracket: end" liquid_p_bracket=(sl_outer.p,sl.p) vapour_p_bracket=(sv.p,sv_outer.p) liquid_v_bracket=(sl_outer.v,sl.v) vapour_v_bracket=(sv.v,sv_outer.v) spinodal_p=(p̄l,p̄v)

    return (sl_outer,sl),(sv,sv_outer),(sml,smv),(p̄l,p̄v)
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
        return xs + xa, evalpoly(xs,poly)
    end
    return zero(hi)/zero(hi),zero(hi)/zero(hi)
end

function coexistence_interpolate_volume(sa,sb,p,liquid)
    if sa.p > sb.p
        slo,shi = sb,sa
    else
        slo,shi = sa,sb
    end
    cp = -slo.p

    dp = shi.p - slo.p
    βp = (p - slo.p)/dp

    if !(zero(βp) <= βp <= one(βp))
        throw(DomainError(lazy"Clapeyron.coexistence_interpolate_volume: pressure $(p) ouside bounds ($(slo.p),$(shi.p))"))
    end

    if !liquid #todo: work in density units
        if βp < 0.8 && abs(shi.dpdv) > 1
            #work in v(p) units and evaluate directly
            poly = Solvers.hermite3_poly(slo.p,shi.p,slo.v,shi.v,1/slo.dpdv,1/shi.dpdv)
            v = evalpoly(p + cp,poly)
        else
            #work in p(v) units and invert via root finding
            poly = Solvers.hermite3_poly(shi.v,slo.v,shi.p,slo.p,shi.dpdv,slo.dpdv)
            k0,k1,k2,k3 = poly
            poly2 = (k0 - p,k1,k2,k3)
            res = Solvers.real_roots3(poly2)
            cv = -shi.v
            lo = shi.v + cv
            hi = slo.v + cv
            in_bounds,vx = Solvers.root_in_bounds(res,lo,hi)
            if in_bounds
                v = vx - cv
            else
                v = shi.v*βp + slo.v*(1 - βp)
            end
        end
    else #work in volume
        if βp > 0.2 && abs(slo.dpdv) > 1
            #work in v(p) units and evaluate directly
            poly = Solvers.hermite3_poly(slo.p,shi.p,slo.v,shi.v,1/slo.dpdv,1/shi.dpdv)
            v = evalpoly(p + cp,poly)
        else
            #work in p(v) units and invert via root finding
            poly = Solvers.hermite3_poly(shi.v,slo.v,shi.p,slo.p,shi.dpdv,slo.dpdv)
            k0,k1,k2,k3 = poly
            poly2 = (k0 - p,k1,k2,k3)
            res = Solvers.real_roots3(poly2)
            cv = -shi.v
            lo = shi.v + cv
            hi = slo.v + cv
            in_bounds,vx = Solvers.root_in_bounds(res,lo,hi)
            if in_bounds
                v = vx - cv
            else
                v = shi.v*βp + slo.v*(1 - βp)
            end
        end
    end
    return v
end