
#=

------------- Coexistence functions --------------------------

Their main use is to bracket a range where there exists a maxwell loop.

=#

function coexistence_volumes_bracketed(model, T, z, _v_lb, _v_ub, vm0 = nothing;
        tol_v    = 1e-12,
        tol_dpdv = 1e-8,
        maxiter  = 100,
        branches = :all)

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

    if isnan(sl.p) || isnan(sv.p)
        return sl,sv,s0,:failure
    end

    # sl is unstable, sv is stable
    if sl.dpdv > 0 && sv.dpdv < 0
        is_vapour(branches) && (return s0,sv,sl,:vapour)
        vl_new = coexistence_expand_volume_bracket(model,T,z,sl.v,sl.dpdv,:liquid)
        sl_new = p_state(vl_new)
        return sl_new,sv,sl,:eq
    end

    # sl is unstable, sv is stable
    if sv.dpdv > 0 && sl.dpdv < 0
        is_liquid(branches) && (return sl,s0,sv,:liquid)
        vv_new = coexistence_expand_volume_bracket(model,T,z,sv.v,sv.dpdv,:vapour)
        sv_new = p_state(vv_new)
        return sl,sv_new,sv,:eq
    end

    # sl is unstable, sv is unstable (??)
    if sv.dpdv > 0 && sl.dpdv > 0
        vv_new = coexistence_expand_volume_bracket(model,T,z,sv.v,sv.dpdv,:vapour)
        vl_new = coexistence_expand_volume_bracket(model,T,z,sl.v,sl.dpdv,:liquid)
        sv_new = p_state(vv_new)
        sl_new = p_state(vl_new)
        smid = ifelse(abs(sl.d2pdv2) < (sv.d2pdv2) ,sl,sv)
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

    for i in 1:maxiter

        lnv_converged,_ = Solvers.convergence(log(sl.v),log(sv.v),zero(TT),1e-12)
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
        if isnan(vm)
            strikes += 1
            vm = sqrt(sl.v*sv.v)
        end

        sm = p_state(vm)

        strikes >= 4 && break

        #return if we found our middle point.
        if sm.dpdv > 0
            status = :eq
            break
        end

        #update bounds
        #Anderson and Prausnitz criteria to distinguish between liquid and gas, good enough here
        #for an stricter criteria, pip was used, but it requires evaluating the EoS again.
        d2pdrho2_m = sm.v * sm.v * sm.v *(sm.v * sm.d2pdv2 + 2 * sm.dpdv) / (∑z*∑z)

        if d2pdrho2_m > zero(TT)
            sl = sm
        else
            sv = sm
        end
    end

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

function coexistence_expand_volume_bracket(model,T,z,v,dpdv,phase)

    #generate a bracket point that actually has dpdv < 0

    TT = Base.promote_eltype(model,T,z,v,dpdv)
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

function coexistence_get_good_bracket(model,T,z,s0l,s0v,s0m)
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

    function p_state(xx)
        _f, _df, _d2f = p∂p∂2p(model, xx, T, z)
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

    #stable outer brackets
    sl_outer,sv_outer = sl,sv

    #estimations of the spinodal pressures
    p̄l,p̄v = TT(NaN),TT(NaN)

    done_l = sl.p < 0
    done_v = false
    cross = false #we stop once we find crossing volumes and the inner bounds cross the middle spinodal pressure.
    #Step 1: find crossing inner bounds, such as sl.p < sv.p. sl also stops if sl.p < 0
    for i in 1:60
        #do at least one iteration to get the spinodal estimation
        cross && done_l && done_v && i > 1 && break
        if !done_l
            poly_l = hermite_from_state(sl,sml) #hermite-3 interpolation to find spindodal
            dpoly_l = Solvers.polyder(poly_l) #dpdv polynomial is a quadratic, the root inside the bounds is the spindodal prediction
            _,vl1,vl2 = Solvers.real_roots2(dpoly_l)
            cl = -sl.v #translation from hermite implementation.
            vsl = ifelse((sl.v + cl) <= vl1 <= (sml.v + cl),vl1,vl2)
            p̄l = evalpoly(vsl,poly_l)
            ssl = p_state(vsl - cl)
            if ssl.dpdv > 0
                vxl1 = sl.v
                vxl2 = sqrt(vl1*ssl.v)
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
            poly_v = hermite_from_state(smv,sv) #hermite-3 interpolation to find spindodal
            dpoly_v = Solvers.polyder(poly_v) #dpdv polynomial is a quadratic, the root inside the bounds is the spindodal prediction
            _,vv1,vv2 = Solvers.real_roots2(dpoly_v)
            cv = -smv.v #translation from hermite implementation.
            vsv = ifelse((smv.v + cv) <= vv1 <= (sv.v + cv),vv1,vv2)
            p̄v = evalpoly(vsv,poly_v)
            ssv = p_state(vsv - cv)
            if ssv.dpdv > 0
                vxv1 = sv.v
                vxv2 = sqrt(ssv.v*vv1)
                smv = ssv
            else
                vxv1 = sv.v
                vxv2 = ssv.v
                sv = ssv
            end
            lnvv_converged,_ = Solvers.convergence(log(vxv1),log(vxv2),zero(TT),1e-12)

            if lnvv_converged
                done_v = true
            end
        end

        if sl.p < sv.p
            cross = true
        end
        
        if cross
            p̄m = 0.5*(max(p̄l,zero(TT)) + p̄v)
            sl.p < p̄m && (done_l = true)
            sv.p > p̄m && (done_v = true)
        end
        @show cross,done_l,done_v,i
    end

    
    #Step 2: find crossing outer bounds, such as sl.p > p̄v, sv.p < p̄l (*)
    # on the case of negative p̄l, we proceed to use the zero volume approximation.

    if sl_outer.p < p̄v
        for i in 1:100
            sl_outer.p > p̄v && break
            vᵢ,pᵢ,dpdvᵢ = sl_outer.v,sl_outer.p,sl_outer.dpdv
            Δᵢ = (1.01p̄v-pᵢ)/(vᵢ*dpdvᵢ) #1.01 is to allow passing p̄v without traps.
            sl_outer = p_state(vᵢ*exp(Δᵢ)) #logVᵢ = logVᵢ + Δᵢ
        end
    end

    RT = Rgas(model)*T
    ∑z = sum(z)
    nRT = ∑z*RT
    if sl_outer.p < 0 || p̄l < 0
        #zero_pressure approximation.
        vl0 = coexistence_interpolate_volume(sl_outer,sl,zero(TT))
        ares_l0 = a_res(model, vl0, T, z)
        lnϕ_l0 = ares_l0 - 1 + log(nRT/vl0)
        p̄l_mod = exp(lnϕ_l0)
    else
        p̄l_mod = p̄l
    end

    if sv_outer.p > p̄l_mod
        #just use ideal gas approximation, after all, with the other bracket we can generate a good EoS in density
        sv_outer = p_state(nRT/p̄l_mod)
    end

    return (sl_outer,sl),(sv,sv_outer),(sml,smv),(p̄l,p̄v)
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