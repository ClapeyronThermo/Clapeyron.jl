normalize_property(model,prop,z,property::F) where F = prop,property
normalize_property(model,prop,z,property::typeof(molar_density)) = sum(z)/prop,volume
normalize_property(model,prop,z,property::typeof(mass_density)) = molecular_weight(model,z)/prop,volume

function x0_edge_temperature(model,p,z,pure = split_pure_model(model))
    dPdTsat = extended_dpdT_temperature.(pure,p)
    T_bubble = antoine_bubble_solve(dPdTsat,p,z)
    T_dew = antoine_dew_solve(dPdTsat,p,z)
    return (T_bubble,T_dew),dPdTsat
end

"""
    edge_temperature(model,p,z,v0 = nothing)

Calculates the temperature at which two fluid phases have the same Gibbs energy and temperature at the specified pressure `p`.

Returns a tuple, containing:
- Edge Temperature `[K]`
- Liquid volume of edge Point `[m³]`
- Vapour volume at edge Point `[m³]`
"""
function edge_temperature(model,p,z,v0 = nothing)
    edge,crit,status = _edge_temperature(model,p,z,v0)
    return edge
end

function edge_temperature(model,p) 
    check_arraysize(model,SA[1.0])
    saturation_temperature(model,p)
end

function try_2ph_edge_temperature2(model,p,z,v10::R,v20::R,T10::R,T20::R,Ts::R) where R
    f = WithContext(edge_temperature_objective2(model,model,p,Ts,Ts,z),∂Tag{∂₁f}())
    V0 = SVector(log(v10),log(v20),T10,T20)
    sol = Solvers.nlsolve2(f,V0,Solvers.Newton2Var())
    v1 = exp(sol[1])
    v2 = exp(sol[2])
    T_eq = 0.5*(sol[3] + sol[4])
    result = (T_eq,v1,v2)
    valid = check_valid_sat_pure(model,p,v1,v2,T_eq,z)
    return result,valid
end

function _edge_temperature(model,p,z,v0 = nothing)
    if v0 == nothing
        vv0,_ = x0_edge_temperature(model,p,z)
    else
        vv0 = (v0[1],v0[2])
    end
    _T1 = vv0[1]
    _T2 = vv0[2]
    Tmin,Tmax = minmax(_T1,_T2)
    n = sum(z)
    v_Tmin = volume(model,p,Tmin,z,phase = :l)
    v_Tmax = volume(model,p,Tmax,z,phase = :v)
    _Ts = 0.5*(_T1 + _T2)
    pp,v10,v20,T10,T20,Ts,_0 = promote(p,v_Tmin,v_Tmax,Tmin,Tmax,_Ts,0.0)
    nan = _0/_0
    fail = (nan,nan,nan)
    _is_positive((v10,v20,T10,T20)) || return fail,fail,:failure

    res0,valid0 = try_2ph_edge_temperature2(model,pp,z,v10,v20,T10,T20,Ts)
    valid0 && (return res0,fail,:success)
    T_eq,vl1,_ = res0

    if isfinite(T_eq)
        res2,_,_ = _edge_pressure(model,T_eq,z,(0.9*p,1.1*p),false)
        p2,v12,v22 = res2
        dpdT = dpdT_saturation(model,model,v12,v22,T_eq,z,z)
        dTinvdlnp = -p2/(dpdT*T_eq*T_eq)
        Δlnp = log(p/p2)
        Tinv0 = 1/T_eq
        Tinv = Tinv0 + dTinvdlnp*Δlnp
        T3 = 1/Tinv
        resx,_,_ = _edge_pressure(model,T3,z,(0.9*p,1.1*p),false)
        _,v1x,v2x = resx
        res1,valid1 = try_2ph_edge_temperature2(model,pp,z,v1x,v2x,T3,T3,Ts)
        valid1 && (return res1,fail,:success)
        if isfinite(res1[2])
            vl1 = res1[2]
        end
    end

    #fail when calculating edge temperature, this happens near the (mechanical) critical point
    if isfinite(T_eq)
        Tr = T_eq/T_scale(model,z)
        vlog = log10(vl1)
        crit = mechanical_critical_point(model,z,(Tr,vlog))
    else
        crit = mechanical_critical_point(model,z)
    end

    Tc,Pc,Vc = crit
    !isfinite(Pc) && return fail,fail,:failure
    Pc <= p && return fail,crit,:supercritical
    T_extrapolated = critical_tsat_extrapolation(model,p,Tc,Pc,Vc*sum(z),z)
    _,vlc,vvc = x0_sat_pure_crit_info(model,T_extrapolated,(Tc,Pc,Vc),z)
    res_crit,valid_crit = try_2ph_edge_temperature2(model,pp,z,vlc,vvc,T_extrapolated,T_extrapolated,Tc)
    valid_crit && (return res_crit,fail,:success)
    return fail,fail,:failure
end


"""
    edge,fa,fb = FindEdge(f::Function,a,b)

Finds approx singularity location in range `a`,`b` for function `f`. There should be only 1 singularity in [`a`,`b`].
Returns the edge point `edge`, and the values at both sides of the edge, sorted such as `a < b`.
"""
function FindEdge(f::T,a,b) where T
    fa,fb = f(a),f(b)
    return FindEdge(f,a,b,fa,fb)
end

function FindEdge(f::T,_a,_b,_fa,_fb) where T
    @assert _a <= _b
    a,b,fa,fb = promote(_a,_b,_fa,_fb)
    for i in 1:40
        isapprox(a,b,rtol=1e-10,atol = 1e-10) && return a,fa,fb
        c = 0.5*(a+b)
        fc = f(c)
        ∇fa,∇fc = (fc - fa)/(c - a),(fb - fc)/(b - a)
        if abs(∇fc) > abs(∇fa)
            a = c
            fa = fc
        else
            b = c
            fb = fc
        end
    end
    nan = zero(a)/zero(a)
    return nan,nan,nan
end

function bubble_temperature_tproperty_method(model,p,T0,z,dPdT)
    y0 = z .* antoine_pressure.(dPdT,T0)
    y0 ./= sum(y0)
    _,T,_,y,vl0,vv0 = improve_bubbledew_suggestion(model,p,T0,z,y0,FugEnum.BUBBLE_TEMPERATURE,FillArrays.Trues(length(z)),false)
    return ChemPotBubbleTemperature((vl0,vv0),T,y,nothing,0.0,1e-8,1e-12,1000,false)
end

function dew_temperature_tproperty_method(model,p,T0,z,dPdT)
    x0 = z ./ antoine_pressure.(dPdT,T0)
    x0 ./= sum(x0)
    _,T,x,_,vl0,vv0 = improve_bubbledew_suggestion(model,p,T0,x0,z,FugEnum.DEW_TEMPERATURE,FillArrays.Trues(length(z)),false)
    return ChemPotDewTemperature((vl0,vv0),T,x,nothing,0.0,1e-8,1e-12,1000,false)
end

@noinline function Xproperty_verbose(type,v1 = nothing,s1 = nothing)
    #error
    type == :error_Tprop    && @error "Tproperty calculation failed"
    type == :error_Pprop    && @error "Pproperty calculation failed"

    #crit
    type == :Pc             && @info "mechanical critical pressure:           $v1"
    type == :Tc             && @info "mechanical critical temperature:        $v1"
    type == :Vc             && @info "mechanical critical molar volume        $v1"
    
    #puresat branch
    type == :satmin_T       && @info "minimum saturation temperature:         $v1"
    type == :satmax_T       && @info "maximum saturation temperature:         $v1"
    type == :satmin_p       && @info "minimum saturation pressure:            $v1"
    type == :satmax_p       && @info "maximum saturation pressure:            $v1"
    type == :satmin_x       && @info "property at minimum sat point:          $v1"
    type == :satmax_x       && @info "property at maximum sat point:          $v1"
    type == :puresat_T      && @info "temperature($s1) outside pure fluid saturation boundaries ($v1)"
    type == :puresat_p      && @info "pressure($s1) outside pure fluid saturation boundaries ($v1)"
    
    #edge
    type == :edge_fail_T    && @warn "failure to calculate edge point, trying to solve using Clapeyron.T_scale(model,z)"
    type == :edge_fail_p    && @warn "failure to calculate edge point, trying to solve using Clapeyron.p_scale(model,z)"
    type == :edge_liq       && @info "property at liquid edge:                $v1"
    type == :edge_vap       && @info "property at vapour edge:                $v1"
    type == :edge_T         && @info "temperature at edge point:              $v1"
    type == :edge_p         && @info "pressure at edge point:                 $v1"

    #inside edge messages
    type == :prop_in_edge   && @info "property between the liquid and vapour edges, in the phase change region"
    type == :prop_may_vap   && @info "property in equilibria, mostly liquid with vapour, checking dew point to improve initial point"
    type == :prop_may_liq   && @info "property in equilibria, mostly vapour with liquid, checking bubble point to improve initial point"

    #outside eq
    type == :outside_eq     && @info "$s1 temperature($v1) outside the phase change region"

    #maybe_bubble, maybe_dew
    type == :dewt           && @info "temperature at dew point:               $v1"
    type == :dewp           && @info "pressure at dew point:                  $v1"
    type == :bubblet        && @info "temperature at bubble point:            $v1"
    type == :bubblep        && @info "pressure at bubble point:               $v1"
    type == :dewx           && @info "property at dew point:                  $v1"
    type == :bubblex        && @info "property at bubble point:               $v1"

    #additional messages for supercritical and pseudo‑phase branches
    type == :bubble_volume  && @info "molar volume at bubble point:           $v1"
    type == :dew_volume     && @info "molar volume at dew point:              $v1"
    type == :volume_at_Tprop&& @info "molar volume at temperature(property):  $v1"
    type == :pseudo_critical&& @info "pseudo-critical temperature($s1) in phase change region (between critical and $v1 points)"
    type == :pseudo_vapour  && @info "pseudo-vapour temperature($s1) in phase change region (between edge and dew point)"
    type == :pseudo_liquid  && @info "pseudo-liquid temperature($s1) in phase change region (between edge and bubble point)"
    type == :pure_over_pc   && @info "pressure is above critical pressure"
    type == :pure_over_tc   && @info "temperature is above critical temperature"
    type == :pure_p_over_sat   && @info "pressure($s1) > saturation pressure"
    type == :pure_p_under_sat  && @info "pressure($s1) < saturation pressure"
    type == :pseudo_critical_p && @info "pseudo-critical pressure($s1) in phase change region (between critical and $v1 points)"
    type == :pseudo_vapour_p   && @info "pseudo-vapour pressure($s1) in phase change region (between edge and dew point)"
    type == :pseudo_liquid_p   && @info "pseudo-liquid pressure($s1) in phase change region (between edge and bubble point)"
end

"""
    Tproperty(model::EoSModel,p,prop,z::AbstractVector,property = enthalpy;rootsolver = Roots.Order0(),phase =:unknown,abstol = 1e-15,reltol = 1e-15, verbose = false)

Given `p` and any other bulk property `prop` calculated via `property`, returns the required temperature `T` such that `property(model,p,T,z,phase) = prop`.

Not all cases of pressure will work as `Clapeyron.bubble_temperature(model,p,z)` and `Clapeyron.dew_temperature(model,p,z)` does not always find a correct starting point.
"""
function Tproperty(model::EoSModel,p,prop,z = SA[1.0],
                    property::TT = enthalpy;
                    rootsolver = Roots.Order0(),
                    phase =:unknown,
                    abstol = 1e-15,
                    reltol = 1e-15,
                    T0 = nothing,
                    verbose = false,
                    threaded = true) where TT
    check_arraysize(model,z)
    cached_model = __tpflash_cache_model(model,p,NaN,z,:vle)
    T,st = _Tproperty(cached_model,p,prop,z,property;rootsolver,phase,abstol,reltol,verbose,threaded,T0)
    return T
end

function __Tproperty_check(res,verbose,Tother = zero(res[1])/zero(res[1]))
    T,st = res
    if st == :failure
        verbose && Xproperty_verbose(:error_Tprop)
        return Tother,st
    end
    return T,st
end

function _Tproperty(model::EoSModel,p,_prop,z = SA[1.0],
                  _property::TT = enthalpy;
                  rootsolver = Roots.Order0(),
                  phase =:unknown,
                  abstol = 1e-15,
                  reltol = 1e-15,
                  T0 = nothing,
                  verbose = false,
                  threaded = true) where TT

    prop,property = normalize_property(model,_prop,z,_property)

    if length(z) == 1
        check_arraysize(model,z)
        zz = SA[z[1]]
        Tnc1,stnc1,_ = Tproperty_pure(fluid_model(model),p,prop,zz,property,rootsolver,phase,abstol,reltol,verbose,threaded,T0)
        return __Tproperty_check((Tnc1,stnc1),verbose)
    end

    if T0 !== nothing
        if is_unknown(phase)
            phase_T0 = identify_phase(model,p,T0,z)
        else
            phase_T0 = phase
        end
        res = __Tproperty(model,p,prop,z,property,rootsolver,phase_T0,abstol,reltol,threaded,T0)
        return __Tproperty_check(res,verbose)
    end

    n = sum(z)
    v0_edge,dpdT = x0_edge_temperature(model,p,z)
    sol_options = (abstol,reltol,rootsolver,verbose)

    if !isfinite(v0_edge[1]) || !isfinite(v0_edge[2]) || !isfinite(prop)
        return __Tproperty_check((NaN*v0_edge[1]*prop,:failure),verbose)
    end

    #check pure saturation envelopes
    Tmin_sat,Tmax_sat = extrema(xx -> T_from_dpdT(xx,p),dpdT)

    update_temperature!(model,Tmin_sat)
    prop_puresat_l = property(model,p,Tmin_sat,z,phase = :l)
    update_temperature!(model,Tmax_sat)
    prop_puresat_v = property(model,p,Tmax_sat,z,phase = :v)
    βpuresat = (prop - prop_puresat_l)/(prop_puresat_v - prop_puresat_l)

    if !(0 <= βpuresat <= 1)  #TODO: check if this is valid
        verbose && Xproperty_verbose(:satmin_T,Tmin_sat)
        verbose && Xproperty_verbose(:satmax_T,Tmax_sat)
        verbose && Xproperty_verbose(:satmin_x,prop_puresat_l)
        verbose && Xproperty_verbose(:satmax_x,prop_puresat_v)

        phase_by_puresat = βpuresat > 1 ? :vapour : :liquid
        if verbose
            Xproperty_verbose(:satmin_T,Tmin_sat)
            Xproperty_verbose(:satmax_T,Tmax_sat)
            Xproperty_verbose(:satmin_x,prop_puresat_l)
            Xproperty_verbose(:satmin_x,prop_puresat_v)
            Xproperty_verbose(:puresat_T,phase_by_puresat,property)
        end
        T_by_puresat = βpuresat > 1 ? Tmax_sat : Tmin_sat
        res_puresat = __Tproperty(model,p,prop,z,property,rootsolver,phase_by_puresat,abstol,reltol,threaded,T_by_puresat)
        return __Tproperty_check(res_puresat,verbose)
    end

    T0_bubble,T0_dew = v0_edge
    edge,crit,status = _edge_temperature(model,p,z,v0_edge)
    T_edge,v_l,v_v = edge

    if is_liquid(phase) || is_vapour(phase)
        if status == :supercritical || status == :failure
            T0_with_phase = is_liquid(phase) ? Tmin_sat : Tmax_sat
        else
            T0_with_phase = is_liquid(phase) ? 0.5*(T_edge + Tmin_sat) : 0.5*(T_edge + Tmax_sat)
        end
        res_with_phase = __Tproperty(model,p,prop,z,property,rootsolver,phase,abstol,reltol,false,T0_with_phase)
        return __Tproperty_check(res_with_phase,verbose)
    end

    if status == :supercritical
        crit_cache = (v0_edge,dpdT,crit)
        return Tproperty_supercritical(model,p,prop,z,property,crit_cache,sol_options)
    end

    if status == :failure
        verbose && Xproperty_verbose(:edge_fail_T)
        res = __Tproperty(model,p,prop,z,property,rootsolver,phase,abstol,reltol,threaded,T_scale(model,z))
        res[2] == :failure && return __Tproperty_check(res,verbose)
        return __Tproperty_check(res,verbose)
    end
    edge_cache = (v0_edge,dpdT,edge)
    T0x,new_phase,prop_edge,success = Tproperty_refine_edge(model,p,prop,z,property,edge_cache,sol_options)
    success && return T0x,new_phase

    if is_vapour(new_phase)
        Tsx,stx,success = Tproperty_maybe_vapour(model,p,prop,z,property,edge_cache,sol_options,prop_edge)
        success && return Tsx,stx
        T0x = Tsx
    end

    if is_liquid(new_phase)
        Tsx,stx,success = Tproperty_maybe_liquid(model,p,prop,z,property,edge_cache,sol_options,prop_edge)
        success && return Tsx,stx
        T0x = Tsx
    end

    verbose && Xproperty_verbose(:outside_eq, new_phase, property)
    res = __Tproperty(model,p,prop,z,property,rootsolver,new_phase,abstol,reltol,threaded,T0x)
    return __Tproperty_check(res,verbose)
end

function Tproperty_supercritical(model,p,prop,z,property,cache,sol_options,phase)
    v0_edge,dpdT,crit = cache
    abstol,reltol,rootsolver,verbose = sol_options
    T0_bubble,T0_dew = v0_edge
    Tc,Pc,Vc = crit
    n = sum(z)

    if is_liquid(phase) || is_vapour(phase)
        res_with_phase = __Tproperty(model,p,prop,z,property,rootsolver,phase,abstol,reltol,false,Tc)
        return __Tproperty_check(res_with_phase,verbose)
    end

    if verbose
        Xproperty_verbose(:Tc,Tc)
        Xproperty_verbose(:Pc,Pc)
        Xproperty_verbose(:Vc,Vc)
    end

    res = __Tproperty(model,p,prop,z,property,rootsolver,:vapour,abstol,reltol,false,Tc)
    res[2] == :failure && return __Tproperty_check(res,verbose)

    #instead of calculating the mixture critical point, we just suppose
    #that all volumes between the bubble and dew volumes evaluated at T = Tc (or P = Pc)
    #are equilibrium ones
    #TODO: we could calculate dvsatdP (or dvsatdT) and estimate a line instead of a vertical threshold
    Tx = res[1]
    Vx = volume(model,p,Tx,z,vol0 = Vc*n)/n

    if Vx <= Vc
        bubble_method_crit = bubble_temperature_tproperty_method(model,Pc,Tc,z,dpdT)
        Tsat,Vsat,_,_ = bubble_temperature(model,Pc,z,bubble_method_crit)
        satpoint = "bubble"
        verbose && Xproperty_verbose(:bubble_volume,Vsat)
    else
        dew_method_crit = dew_temperature_tproperty_method(model,Pc,Tc,z,dpdT)
        Tsat,_,Vsat,_ = dew_temperature(model,Pc,z,dew_method_crit)
        satpoint = "dew"
        verbose && Xproperty_verbose(:dew_volume,Vsat)
    end
    verbose && Xproperty_verbose(:volume_at_Tprop,Vx)

    βx = (Vx - Vsat)/(Vc - Vsat)
    if 0 <= βx <= 1
        verbose && Xproperty_verbose(:pseudo_critical, satpoint, property)
        return Tx,:eq
    end

    verbose && @info "temperature(property) in the critical pseudo-$(string(res[2])) branch, outside the phase change region"
    return res
end

function Tproperty_refine_edge(model,p,prop,z,property,cache,sol_options)
    v0_edge,dpdT,edge = cache
    abstol,reltol,rootsolver,verbose = sol_options
    T0_bubble,T0_dew = v0_edge
    T_edge,v_l,v_v = edge

    if has_a_res(model)
        prop_l = spec_to_vt(model,v_l,T_edge,z,property)
        prop_v = spec_to_vt(model,v_v,T_edge,z,property)
    else
        prop_l = property(model,p,T_edge,z,phase = :l)
        prop_v = property(model,p,T_edge,z,phase = :v,vol0 = v_v)
    end

    if verbose
        Xproperty_verbose(:edge_liq,prop_l)
        Xproperty_verbose(:edge_vap,prop_v)
        Xproperty_verbose(:edge_T,T_edge)
    end

    β_edge = (prop - prop_l)/(prop_v - prop_l)

    #we are inside equilibria.
    if 0 <= β_edge <= 1
        verbose && Xproperty_verbose(:prop_in_edge)
        T_edge_interp = β_edge*T0_dew + (1 - β_edge)*T0_bubble
        β_T_edge = (T_edge - T0_dew)/(T0_bubble - T0_dew)
        
        ϕ = 0.3 #P_edge is in the center of the bubble and dew approximations, return P_edge
        if ϕ <= β_T_edge <= (1 - ϕ)
            return T_edge_interp,:eq,one(β_edge)*prop,true
        end
            
        if β_T_edge < 0.3
            T0x = T0_bubble
            verbose && Xproperty_verbose(:prop_may_vap)
            #we search between the liquid edge and the dew temperature
            update_temperature!(model,T0x)
            prop_edge,new_phase = property(model,p,T0_bubble,z,phase = :l),:vapour
        else
            T0x = T0_dew
            verbose && Xproperty_verbose(:prop_may_liq)
            #we search between the vapour edge and the bubble temperature
            update_temperature!(model,T0x)
            prop_edge,new_phase = property(model,p,T0_dew,z,phase = :v),:liquid
        end
    else
        T0x = T_edge
        if β_edge > 1
            prop_edge,new_phase = prop_v,:vapour
        else
            prop_edge,new_phase = prop_l,:liquid
        end
    end

    return T0x,new_phase,prop_edge,false
end

function Tproperty_maybe_vapour(model,p,prop,z,property,cache,sol_options,prop_edge)
    v0_edge,dpdT,edge = cache
    abstol,reltol,rootsolver,verbose = sol_options
    T0_bubble,T0_dew = v0_edge
    T_edge,_,_ = edge
    n = sum(z)

    dew_method = dew_temperature_tproperty_method(model,p,T0_dew,z,dpdT)
    dew = dew_temperature(model,p,z,dew_method)
    T_dew,_,v_dew,_ = dew
    update_temperature!(model,T_dew)
    if has_a_res(model)
        prob_dew = spec_to_vt(model,v_dew*n,T_dew,z,property)
    else
        prob_dew = property(model,p,T_dew,z,phase = :vapour,vol0 = n*v_dew)
    end

    verbose && Xproperty_verbose(:dewt,T_dew)
    verbose && Xproperty_verbose(:dewx,prob_dew)
    β_dew = (prop - prop_edge)/(prob_dew - prop_edge)
    T_interp = β_dew*T_dew + (1 - β_dew)*T_edge

    if 0 < β_dew < 1
        verbose && Xproperty_verbose(:pseudo_vapour, property)
        T1,st1 = __Tproperty_check((T_interp,:eq),verbose,T_edge)
        return T1,st1,true
    else
        return T_dew,:vapour,false
    end     
end

function Tproperty_maybe_liquid(model,p,prop,z,property,cache,sol_options,prop_edge)
    v0_edge,dpdT,edge = cache
    abstol,reltol,rootsolver,verbose = sol_options
    T0_bubble,T0_dew = v0_edge
    T_edge,_,_ = edge
    n = sum(z)

    bubble_method = bubble_temperature_tproperty_method(model,p,T0_bubble,z,dpdT)    
    bubble = bubble_temperature(model,p,z,bubble_method)
    T_bubble,v_bubble,_,_ = bubble
    update_temperature!(model,T_bubble)
    if has_a_res(model)
        prob_bubble = spec_to_vt(model,v_bubble*n,T_bubble,z,property)
    else
        prob_bubble = property(model,p,T_bubble,z,phase = :l)
    end

    verbose && Xproperty_verbose(:bubblet,T_bubble)
    verbose && Xproperty_verbose(:bubblex,prob_bubble)

    β_bubble = (prop - prop_edge)/(prob_bubble - prop_edge)
    T_interp = β_bubble*T_bubble + (1 - β_bubble)*T_edge

    if 0 < β_bubble < 1
        verbose && Xproperty_verbose(:pseudo_liquid, property)
        T1,st1 = __Tproperty_check((T_interp,:eq),verbose,T_edge)
        return T1,st1,true
    else
        return T_bubble,:liquid,false
    end   
end

Tproperty_pure(model,p,x,z,property::F,phase,T0,verbose) where F = Tproperty_pure(model,p,x,z,property,Roots.Order0(),phase,1e-15,1e-15,verbose,false,T0)

function Tproperty_pure(model,p,x,z,property::F,rootsolver,phase,abstol,reltol,verbose,threaded,T0) where F
    TT = Base.promote_eltype(model,p,x,z)
    nan = zero(TT)/zero(TT)
    ∑z = sum(z)
    x1 = SVector(1.0*one(∑z))

    if !isnothing(T0)
        if is_unknown(phase)
            phase_T0 = identify_phase(model,p,T0,z)
        else
            phase_T0 = phase
        end
        T_init,_ = __Tproperty(model,p,x,z,property,rootsolver,phase_T0,abstol,reltol,false,T0)
        return T_init,phase_T0,(nan,nan,nan)
    end

    sat,crit,status = _extended_saturation_temperature(model,p)

    if status == :failure
        verbose && Xproperty_verbose(:error_Tprop)
        return nan,:failure,(nan,nan,nan)
    end

    if status == :supercritical
        verbose && Xproperty_verbose(:pure_over_Pc)
        Tc,Pc,Vc = crit      
        Tcrit0 = TT(1.001Tc) #some eos have problems at exactly the critical point (SingleFluid("R123"))
        Tsc,st_sc = __Tproperty(model,p,x,z,property,rootsolver,:liquid,abstol,reltol,false,Tcrit0)
        return Tsc,st_sc,(nan,nan,nan)
    end

    Ts,vl,vv = TT.(sat)
    
    xl = ∑z*spec_to_vt(model,vl,Ts,x1,property)
    xv = ∑z*spec_to_vt(model,vv,Ts,x1,property)
    βv = (x - xl)/(xv - xl)

    if !isfinite(βv)
        verbose && Xproperty_verbose(:error_Tprop)
        return nan,:failure,nan
    elseif βv < 0 || βv > 1
        phase0 = βv < 0 ? :liquid : :vapour
        is_liquid(phase0) && verbose && @info "temperature($property) < saturation temperature"
        is_vapour(phase0) && verbose && @info "temperature($property) > saturation temperature"
        T1ph,st1ph = __Tproperty(model,p,x,z,property,rootsolver,phase0,abstol,reltol,threaded,Ts)
        return T1ph,st1ph,(nan,nan,nan)
    else
        verbose && Xproperty_verbose(:prop_in_edge)
        return Ts,:eq,(βv,vl,vv)
    end
end

function __Tproperty(model,p,prop,z,property::F,rootsolver,phase,abstol,reltol,threaded,T0) where F
    if is_unknown(phase)
        new_phase = identify_phase(model,p,T0,z)
        if is_unknown(new_phase) #something really bad happened
            _0 = zero(Base.promote_eltype(model,p,prop,z))
            nan = _0/_0
            return nan,:unknown
        end
    else
        new_phase = phase
    end

    _threaded = length(z) == 1 ? false : threaded
    _1 = oneunit(typeof(prop))
    function f(t,tup)
        _prop,_model,_p,_z = tup
        update_temperature!(model,t)
        _1*property(_model,_p,t,_z,phase = new_phase,threaded = _threaded) - _prop
    end

    prob_params = (prop,model,p,z)

    prob = Roots.ZeroProblem(f,_1*T0)
    T = roots_solve_ad(prob,rootsolver,prob_params,atol = abstol,rtol = reltol)
    return T,phase
end

__Tproperty(model,p,prop,z,property::F,phase,T0) where F = __Tproperty(model,p,prop,z,property,Roots.Order0(),phase,1e-15,1e-15,true,T0)
__Tproperty(model,p,prop,z,property::F,phase,T0,verbose::Bool) where F = __Tproperty(model,p,prop,z,property,Roots.Order0(),phase,1e-15,1e-15,verbose,T0)

function __Tproperty(model,p,prop,property::F,rootsolver,phase,abstol,reltol,threaded,T0) where F
    __Tproperty(model,p,prop,SA[1.0],property,rootsolver,phase,abstol,reltol,threaded,T0)
end

# model = PCSAFT(["propane","dodecane"])
# p = 101325; T = 300;
# z = [0.5,0.5]
# h_ = enthalpy(model,p,T,z)
# s_ = entropy(model,p,T,z)
# cp_ = isobaric_heat_capacity(model,p,T,z)
# ρ_ = mass_density(model,p,T,z)
# ic_ = Clapeyron.isentropic_compressibility(model,p,T,z)

# sol1 = Tproperty(model,p,h_,z,enthalpy)
# sol2 = Tproperty(model,p,s_,z,entropy)
# sol3 = Tproperty(model,p,ρ_,z,mass_density)
# sol4 = Tproperty(model,p,ic_,z,isentropic_compressibility)

export Tproperty, edge_temperature