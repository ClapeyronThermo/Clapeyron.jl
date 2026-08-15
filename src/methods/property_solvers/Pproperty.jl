function x0_edge_pressure(model,T,z,pure = split_pure_model(model))
    sat = extended_saturation_pressure.(pure,T)
    n = sum(z)
    p_bubble = sum(z[i]*first(sat[i]) for i in 1:length(model))/n
    p_dew = n/sum(z[i]/first(sat[i]) for i in 1:length(model))
    return (p_bubble,p_dew),sat
end

"""
    edge_pressure(model,T,z,v0 = nothing)

Calculates the pressure at which two fluid phases have the same Gibbs energy and pressure at the specified temperature `T`.

Returns a tuple, containing:
- Edge Pressure `[Pa]`
- Liquid volume of edge Point `[m³]`
- Vapour volume at edge Point `[m³]`
"""
function edge_pressure(model,T,z,v0 = nothing;crit_retry = true)
    edge,crit,status = _edge_pressure(model,T,z,v0,crit_retry)
    return edge
end

function _edge_pressure(model,T,z,v0 = nothing,crit_retry = true)
    if v0 == nothing
        vv0,_ = x0_edge_pressure(model,T,z)
    else
        vv0 = (v0[1],v0[2])
    end
    p1 = vv0[1]
    p2 = vv0[2]
    pmin,pmax = minmax(p1,p2)

    if pmax/pmin > 10
        p_near0,_,_ = x0_sat_pure_near0(model,T,z)
        pmin = 0.9*p_near0
        pmax = 10*p_near0
    end

    v_pmin = volume(model,pmin,T,z,phase = :v)
    v_pmax = volume(model,pmax,T,z,phase = :l)
    TT = T*one(Base.promote_eltype(model,v_pmin,v_pmax,T))
    V01,V02,_0  = promote(v_pmax,v_pmin,zero(TT))
    nan = _0/_0
    fail = (nan,nan,nan)
    _is_positive((v_pmin,v_pmax,T)) || return fail,fail,:failure
    ps,mus = equilibria_scale(model,z)
    edge,valid0 = try_2ph_edge_pressure(model,model,T,V01,V02,ps,mus,z,nothing)
    p_eq,v1,v2 = edge
    valid0 && return edge,fail,:success
    !crit_retry && return fail,fail,:failure
    #fail when calculating edge pressure, this happens near the (mechanical) critical point
    Tr = T/T_scale(model,z)
    vlog = log10(v1)
    crit = mechanical_critical_point(model,z,(Tr,vlog)) #mechanical critical point
    Tc,Pc,Vc = crit
    !isfinite(Tc) && return fail,fail,:failure
    Tc <= T && return fail,crit,:supercritical
    vlc,vvc = critical_vsat_extrapolation(model,T,Tc,Vc,z)
    edge_crit,valid_crit = try_2ph_edge_pressure(model,model,T,vlc,vvc,ps,mus,z,nothing)
    valid_crit && return edge_crit,crit,:success
    return fail,fail,:failure
end

function bubble_pressure_pproperty_method(model,p0,T,z,sat)
    y0 = z .* first.(sat)
    y0 ./= sum(y0)
    p,_,_,y,vl0,vv0 = improve_bubbledew_suggestion(model,p0,T,z,y0,FugEnum.BUBBLE_PRESSURE,FillArrays.Trues(length(z)),false)
    return ChemPotBubblePressure((vl0,vv0),p,y,nothing,0.0,1e-8,1e-12,1000,false)
end

function dew_pressure_pproperty_method(model,p0,T,z,sat)
    x0 = z ./ first.(sat)
    x0 ./= sum(x0)
    p,_,x,_,vl0,vv0 = improve_bubbledew_suggestion(model,p0,T,x0,z,FugEnum.DEW_PRESSURE,FillArrays.Trues(length(z)),false)
    return ChemPotDewPressure((vl0,vv0),p,x,nothing,0.0,1e-8,1e-12,1000,false)
end

"""
    Pproperty(model::EoSModel,T,prop,z = SA[1.0],property::TT = enthalpy;rootsolver = Roots.Order0(),phase =:unknown,abstol = 1e-15,reltol = 1e-15, verbose = false)

Given `T` and any other bulk property `prop` calculated via `property`, returns the required pressure `p` such that `property(model,p,T,z,phase) = prop`.

Not all cases of temperature will work as `Clapeyron.bubble_pressure(model,T,z)` and `Clapeyron.dew_pressure(model,T,z)` does not always find a correct starting point.
"""
function Pproperty(model::EoSModel,T,prop,z = SA[1.0],
                  property::TT = enthalpy;
                  rootsolver = Roots.Order0(),
                  phase =:unknown,
                  abstol = 1e-15,
                  reltol = 1e-15,
                  p0 = nothing,
                  verbose = false,
                  threaded = true) where TT

    cached_model = __tpflash_cache_model(model,NaN,T,z,:vle)
    p,st = _Pproperty(cached_model,T,prop,z,property;rootsolver,phase,abstol,reltol,p0,verbose,threaded)
    return p
end

function __Pproperty_check(res,verbose,p_other = zero(res[1])/zero(res[1]))
    p,st = res
    if st == :failure
        verbose && Xproperty_verbose(:error_Pprop)
        return p_other,st
    end
    return p,st
end

function _Pproperty(model::EoSModel,T,_prop,z = SA[1.0],
                    _property::TT = enthalpy;
                    rootsolver = Roots.Order0(),
                    phase =:unknown,
                    abstol = 1e-15,
                    reltol = 1e-15,
                    p0 = nothing,
                    verbose = false,
                    threaded = true) where TT

    prop,property = normalize_property(model,_prop,z,_property)

    if length(z) == 1
        check_arraysize(model,z)
        zz = SA[z[1]]
        pnc1,ptnc1,_ = Pproperty_pure(fluid_model(model),T,prop,zz,property,rootsolver,phase,abstol,reltol,verbose,threaded,p0)
        return __Pproperty_check((pnc1,ptnc1),verbose)
    end

    if p0 !== nothing
        res = __Pproperty(model,T,prop,z,property,rootsolver,phase_p0,abstol,reltol,threaded,p0)
        return __Pproperty_check(res,verbose)
    end

    n = sum(z)
    v0_edge,pure_sats = x0_edge_pressure(model,T,z)
    sol_options = (abstol,reltol,rootsolver,verbose)

    if !isfinite(v0_edge[1]) || !isfinite(v0_edge[2]) || !isfinite(prop)
        return __Pproperty_check((NaN*v0_edge[1]*prop,:failure),verbose)
    end

    #check pure saturation envelopes
    p_puresat,st_puresat,pb = Pproperty_puresat(model,T,prop,z,property,(v0_edge,pure_sats),sol_options,phase)
    st_puresat != :failure && return (p_puresat,st_puresat)
    pmin_sat,pmax_sat = pb

    #check isogibbs condition ("edge")
    p0_bubble,p0_dew = v0_edge
    edge,crit,status = _edge_pressure(model,T,z,v0_edge)
    P_edge,v_l,v_v = edge
    edge_cache = (v0_edge,pure_sats,edge)

    if is_liquid(phase) || is_vapour(phase)
        if status == :supercritical || status == :failure
            p0_with_phase = is_liquid(phase) ? pmax_sat : pmin_sat
        else
            p0_with_phase = is_liquid(phase) ? 0.5*(P_edge + pmax_sat) : 0.5*(P_edge + pmin_sat)
        end
        res_with_phase = __Pproperty(model,T,prop,z,property,rootsolver,phase,abstol,reltol,false,p0_with_phase)
        return __Pproperty_check(res_with_phase,verbose)
    end

    if status == :supercritical
        crit_cache = (v0_edge,pure_sats,crit)
        return Pproperty_supercritical(model,T,prop,z,property,crit_cache,sol_options,phase)
    end

    if status == :failure
        verbose && Xproperty_verbose(:edge_fail_p)
        res = __Pproperty(model,T,prop,z,property,rootsolver,phase,abstol,reltol,threaded,p_scale(model,z))
        res[2] == :failure && return __Pproperty_check(res,verbose)
        return __Pproperty_check(res,verbose)
    end
    
    p0x,new_phase,prop_edge,success = Pproperty_refine_edge(model,T,prop,z,property,edge_cache,sol_options)
    success && return p0x,new_phase

    if is_vapour(new_phase)
        psx,stx,success = Pproperty_maybe_vapour(model,T,prop,z,property,edge_cache,sol_options,prop_edge)
        success && return psx,stx
        p0x = psx
    end

    if is_liquid(new_phase)
        psx,stx,success = Pproperty_maybe_liquid(model,T,prop,z,property,edge_cache,sol_options,prop_edge)
        success && return psx,stx
        p0x = psx
    end

    verbose && Xproperty_verbose(:outside_eq, new_phase, property)
    res = __Pproperty(model,T,prop,z,property,rootsolver,new_phase,abstol,reltol,threaded,p0x)
    return __Pproperty_check(res,verbose)
end

#we check if the property lies outside the extended bound defined by the extrema of saturation pressures.
function Pproperty_puresat(model,T,prop,z,property,cache,sol_options,phase)
    v0_edge,pure_sats = cache
    abstol,reltol,rootsolver,verbose = sol_options
    pmin_sat,pmax_sat = extrema(first,pure_sats)
    prop_puresat_l = property(model,pmax_sat,T,z,phase = :l)
    prop_puresat_v = property(model,pmin_sat,T,z,phase = :v)
    βpuresat = (prop - prop_puresat_l)/(prop_puresat_v - prop_puresat_l)

    if !(0 <= βpuresat <= 1)  #TODO: check if this is valid
        phase_puresat = βpuresat > 1 ? :vapour : :liquid
        if verbose
            verbose && Xproperty_verbose(:satmin_p,pmin_sat)
            verbose && Xproperty_verbose(:satmax_p,pmax_sat)
            verbose && Xproperty_verbose(:satmin_x,prop_puresat_l)
            verbose && Xproperty_verbose(:satmax_x,prop_puresat_v)
            verbose && Xproperty_verbose(:puresat_p, phase_puresat, property)
        end

        #specified phase is not equal to estimated phase, bail out.
        if is_liquid(phase) && is_vapour(phase_puresat)
            pmin_sat,:failure,(pmin_sat,pmax_sat)
        end

        if is_vapour(phase) && is_liquid(phase_puresat)
            pmin_sat,:failure,(pmin_sat,pmax_sat)
        end

        p_puresat0 = βpuresat > 1 ? pmin_sat : pmax_sat
        res_puresat = __Pproperty(model,T,prop,z,property,rootsolver,phase_puresat,abstol,reltol,threaded,p_puresat0)
        p_puresat,st_puresat = __Pproperty_check(res_puresat,verbose)
        return p_puresat,st_puresat,(pmin_sat,pmax_sat)
    end
    return pmin_sat,:failure,(pmin_sat,pmax_sat)
end

function Pproperty_supercritical(model,T,prop,z,property,cache,sol_options,phase)
    v0_edge,pure_sats,crit = cache
    abstol,reltol,rootsolver,verbose = sol_options
    p0_bubble,p0_dew = v0_edge
    Tc,Pc,Vc = crit
    n = sum(z)

    if verbose
        Xproperty_verbose(:Tc,Tc)
        Xproperty_verbose(:Pc,Pc)
        Xproperty_verbose(:Vc,Vc)
    end

    res = __Pproperty(model,T,prop,z,property,rootsolver,:vapour,abstol,reltol,false,Pc)
    res[2] == :failure && return __Pproperty_check(res,verbose)

    #instead of calculating the mixture critical point, we just suppose
    #that all volumes between the bubble and dew volumes evaluated at T = Tc (or P = Pc)
    #are equilibrium ones
    #TODO: we could calculate dvsatdP (or dvsatdT) and estimate a line instead of a vertical threshold
    px = res[1]
    Vx = volume(model,px,T,z,vol0 = Vc*n)/n

    if Vx <= Vc
        bubble_method_crit = bubble_pressure_pproperty_method(model,Pc,Tc,z,pure_sats)
        Psat,Vsat,_,_ = bubble_pressure(model,Tc,z,bubble_method_crit)
        verbose && Xproperty_verbose(:bubble_volume,Vsat)
    else
        dew_method_crit = dew_pressure_pproperty_method(model,Pc,Tc,z,pure_sats)
        Psat,_,Vsat,_ = dew_pressure(model,Tc,z,dew_method_crit)
        verbose && Xproperty_verbose(:dew_volume,Vsat)
    end

    verbose && Xproperty_verbose(:volume_at_Tprop,Vx)

    βx = (Vx - Vsat)/(Vc - Vsat)
    if 0 <= βx <= 1
        verbose && Xproperty_verbose(:pseudo_critical_p, satpoint, property)
        return px,:eq
    end

    if verbose
        if Vx <= Vc
            Xproperty_verbose(:outside_eq_p, :bubble, property)
        elseif Vx > Vc
            Xproperty_verbose(:outside_eq_p, :dew, property)
        end
    end

    return res
end

function Pproperty_refine_edge(model,T,prop,z,property,cache,sol_options)
    v0_edge,pure_sats,edge = cache
    abstol,reltol,rootsolver,verbose = sol_options
    p0_bubble,p0_dew = v0_edge
    P_edge,v_l,v_v = edge

    prop_l = spec_to_vt(model,v_l,T,z,property)
    prop_v = spec_to_vt(model,v_v,T,z,property)

    verbose && Xproperty_verbose(:edge_liq,prop_l)
    verbose && Xproperty_verbose(:edge_vap,prop_v)
    verbose && Xproperty_verbose(:edge_p,P_edge)

    β_edge = (prop - prop_l)/(prop_v - prop_l)

    #we are inside equilibria.
    if 0 <= β_edge <= 1
        verbose && Xproperty_verbose(:prop_in_edge)
        P_edge_interp = β_edge*p0_dew + (1 - β_edge)*p0_bubble
        β_P_edge = (P_edge - p0_dew)/(p0_bubble - p0_dew)

        ϕ = 0.3 #P_edge is in the center of the bubble and dew approximations, return P_edge
        if ϕ <= β_P_edge <= (1 - ϕ)
            return P_edge_interp,:eq,one(β_edge)*prop,true
        end

        if β_P_edge < 0.3
            p0x = p0_bubble
            verbose && Xproperty_verbose(:prop_may_vap)
            #we search between the liquid edge and the dew pressure
            prop_edge,new_phase = property(model,p0_bubble,T,z,phase = :l),:vapour
        else
            p0x = p0_dew
            verbose && Xproperty_verbose(:prop_may_liq)
            #we search between the vapour edge and the bubble pressure
            prop_edge,new_phase = property(model,p0_dew,T,z,phase = :v),:liquid
        end
    else
        p0x = P_edge
        if β_edge > 1
            prop_edge,new_phase = prop_v,:vapour
        else
            prop_edge,new_phase = prop_l,:liquid
        end
    end

    return p0x,new_phase,prop_edge,false
end

function Pproperty_maybe_vapour(model,T,prop,z,property,cache,sol_options,prop_edge)
    v0_edge,pure_sats,edge = cache
    abstol,reltol,rootsolver,verbose = sol_options
    p0_bubble,p0_dew = v0_edge
    P_edge,_,_ = edge
    n = sum(z)

    dew_method = dew_pressure_pproperty_method(model,p0_dew,T,z,pure_sats)
    dew = dew_pressure(model,T,z,dew_method)
    p_dew,_,v_dew,_ = dew
    prob_dew = spec_to_vt(model,v_dew*n,T,z,property)

    verbose && Xproperty_verbose(:dewp,p_dew)
    verbose && Xproperty_verbose(:dewx,prob_dew)

    β_dew = (prop - prop_edge)/(prob_dew - prop_edge)
    p_interp = exp(β_dew*log(p_dew) + (1 - β_dew)*log(P_edge))

    if 0 < β_dew < 1
        verbose && Xproperty_verbose(:pseudo_vapour_p, property)
        p1,st1 = __Pproperty_check((p_interp,:eq),verbose,P_edge)
        return p1,st1,true
    else
        return p_dew,:vapour,false
    end
end

function Pproperty_maybe_liquid(model,T,prop,z,property,cache,sol_options,prop_edge)
    v0_edge,pure_sats,edge = cache
    abstol,reltol,rootsolver,verbose = sol_options
    p0_bubble,p0_dew = v0_edge
    P_edge,_,_ = edge
    n = sum(z)

    bubble_method = bubble_pressure_pproperty_method(model,p0_bubble,T,z,pure_sats)
    bubble = bubble_pressure(model,T,z,bubble_method)
    p_bubble,v_bubble,_,_ = bubble
    prob_bubble = spec_to_vt(model,v_bubble*n,T,z,property)

    verbose && Xproperty_verbose(:bubblep,p_bubble)
    verbose && Xproperty_verbose(:bubblex,prob_bubble)

    β_bubble = (prop - prop_edge)/(prob_bubble - prop_edge)
    p_interp = exp(β_bubble*log(p_bubble) + (1 - β_bubble)*log(P_edge))

    if 0 < β_bubble < 1
        verbose && Xproperty_verbose(:pseudo_liquid_p, property)
        p1,st1 = __Pproperty_check((p_interp,:eq),verbose,P_edge)
        return p1,st1,true
    else
        return p_bubble,:liquid,false
    end
end

Pproperty_pure(model,T,x,z,property::F,phase,p0,verbose) where F = Tproperty_pure(model,T,x,z,property,Roots.Order0(),phase,1e-15,1e-15,verbose,false,p0)

function Pproperty_pure(model,T,x,z,property::F,rootsolver,phase,abstol,reltol,verbose,threaded,p0) where F
    TT = Base.promote_eltype(model,T,x,z)
    nan = zero(TT)/zero(TT)
    ∑z = sum(z)
    x1 = SA[1.0*one(∑z)]

    if !isnothing(p0)
        p_init,phase_p0 = __Pproperty(model,T,x,z,property,rootsolver,phase,abstol,reltol,false,p0)
        return p_init,phase_p0,(nan,nan,nan)
    end

    sat,crit,status = _extended_saturation_pressure(model,T)

    if status == :failure
        verbose && Xproperty_verbose(:error_Pprop)
        return nan,status,(nan,nan,nan)
    end

    if status == :supercritical
        verbose && Xproperty_verbose(:pure_over_Tc)
        Tc,Pc,Vc = crit
        if p0 !== nothing
            Pcrit0 = TT(p0)
        else
            Pcrit0 = TT(1.001Pc) #some eos have problems at exactly the critical point (SingleFluid("R123"))
        end
        psc,st_sc = __Pproperty(model,T,x,z,property,rootsolver,:liquid,abstol,reltol,false,Pcrit0)
        return psc,st_sc,(nan,nan,nan)
    end

    ps,vl,vv = TT.(sat)

    xl = ∑z*spec_to_vt(model,vl,T,x1,property)
    xv = ∑z*spec_to_vt(model,vv,T,x1,property)
    βv = (x - xl)/(xv - xl)

    if !isfinite(βv)
        verbose && Xproperty_verbose(:error_Pprop)
        return nan,:failure,(nan,nan,nan)
    elseif βv < 0 || βv > 1
        phase0 = βv < 0 ? :liquid : :vapour
        is_liquid(phase0) && verbose && Xproperty_verbose(:pure_p_over_sat, property)
        is_vapour(phase0) && verbose && Xproperty_verbose(:pure_p_under_sat, property)
        p1ph,st1ph = __Pproperty(model,T,x,z,property,rootsolver,phase0,abstol,reltol,threaded,ps)
        return p1ph,st1ph,(nan,nan,nan)
    else
        verbose && Xproperty_verbose(:prop_in_edge)
        return ps,:eq,(βv,vl,vv)
    end
end

function __Pproperty(model,T,prop,z,property::F,rootsolver,phase,abstol,reltol,threaded,p0) where F
    new_phase = is_unknown(phase) ? identify_phase(model,p0,T,z) : phase

    if is_unknown(new_phase) #something really bad happened
        _0 = zero(Base.promote_eltype(model,T,prop,z))
        nan = _0/_0
        return nan,:unknown
    end

    # special case: if property is volume, we can solve pressure directly without root finding
    if property == volume
        px = pressure(model,prop,T,z)
        phasex = VT_identify_phase(model,prop,T,z)
        if is_unknown(phasex)
            return px,:failure
        else
            return px,phasex
        end
    end

    _threaded = length(z) == 1 ? false : threaded
    _1 = oneunit(typeof(prop))
    function f(lnp,tup)
        _prop,_model,_T,_z = tup
        _1*property(_model,exp(lnp),_T,_z,phase = new_phase,threaded = _threaded) - _prop
    end

    prob_params = (prop,model,T,z)

    prob = Roots.ZeroProblem(f, _1*log(p0))
    logp = roots_solve_ad(prob,rootsolver,prob_params,atol = abstol,rtol = reltol)
    return exp(logp),new_phase
end

#=
model = PCSAFT(["propane","dodecane"])
p = 3*101325; T = 300;
z = [0.5,0.5]
h_ = enthalpy(model,p,T,z)
s_ = entropy(model,p,T,z)
ρ_ = mass_density(model,p,T,z)


sol1 = Pproperty(model,T,h_,z,enthalpy)
sol2 = Pproperty(model,T,s_,z,entropy)
sol3 = Pproperty(model,T,ρ_,z,mass_density) =#

export Pproperty, edge_pressure