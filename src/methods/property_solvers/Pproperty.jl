function x0_edge_pressure(model,T,z,pure = split_pure_model(model))
  sat = extended_saturation_pressure.(pure,T)
  n = sum(z)
  p_bubble = sum(z[i]*first(sat[i]) for i in 1:length(model))/n
  p_dew = n/sum(z[i]/first(sat[i]) for i in 1:length(model))
  return (p_bubble,p_dew),sat
end

"""
    edge_pressure(model,p,z,v0 = nothing)

Calculates the pressure at which two fluid phases have the same gibbs and pressure at the specified temperature.

Returns a tuple, containing:
- Edge Pressure `[Pa]`
- Liquid volume of edge Point `[m³]`
- Vapour volume at edge Point `[m³]`
"""
function edge_pressure(model,T,z,v0 = nothing)
  edge,crit,status = _edge_pressure(model,T,z,v0)
  return edge
end

function _edge_pressure(model,T,z,v0 = nothing)
  if v0 == nothing
    vv0,_ = x0_edge_pressure(model,T,z)
  else
    vv0 = (v0[1],v0[2])
  end
  p1 = vv0[1]
  p2 = vv0[2]
  pmin,pmax = minmax(p1,p2)
  v_pmin = volume(model,pmin,T,z,phase = :v)
  v_pmax = volume(model,pmax,T,z,phase = :l)
  f(x) = μp_equality1_p(model,exp(x[1]),exp(x[2]),T,z)
  TT = T*one(Base.promote_eltype(model,v_pmin,v_pmax,T))
  V0 = svec2(log(v_pmax),log(v_pmin),TT)

  _0 = zero(V0[1])
  nan = _0/_0
  fail = (nan,nan,nan)

  _is_positive((v_pmin,v_pmax,T)) || return fail,fail,:failure

  sol = Solvers.nlsolve2(f,V0,Solvers.Newton2Var())
  v1 = exp(sol[1])
  v2 = exp(sol[2])
  p_eq = pressure(model,v2,T,z)
  edge = (p_eq,v1,v2)
  check_valid_sat_pure(model,p_eq,v1,v2,T,z) && return edge,fail,:success

  #fail when calculating edge pressure, this happens near the (mechanical) critical point
  Tr = T/T_scale(model,z)
  vlog = log10(v1)
  crit = mechanical_critical_point(model,z,(Tr,vlog)) #mechanical critical point
  Tc,Pc,Vc = crit

  !isfinite(Tc) && return fail,fail,:failure
  Tc <= T && return fail,crit,:supercritical

  vlc,vvc = critical_vsat_extrapolation(model,T,Tc,Vc,z)
  V1 = SVector(promote(log(vlc),log(vvc)))
  sol1 = Solvers.nlsolve2(f,V1,Solvers.Newton2Var())
  v3 = exp(sol1[1])
  v4 = exp(sol1[2])
  p_eq2 = pressure(model,v2,T,z)
  edge2 = (v3,v4,p_eq2)
  check_valid_sat_pure(model,p_eq2,v3,v4,T,z) && return (edge2,crit,:success)

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
  p,st = _Pproperty(model,T,prop,z,property;rootsolver,phase,abstol,reltol,p0,verbose,threaded)
  return p
end

function __Pproperty_check(res,verbose,p_other = zero(res[1])/zero(res[1]))
  p,st = res
  if verbose && st == :failure
    @error "PProperty calculation failed"
    return p_other,st
  end
  return p,st
end

function _Pproperty(model::EoSModel,T,prop,z = SA[1.0],
          property::TT = enthalpy;
          rootsolver = Roots.Order0(),
          phase =:unknown,
          abstol = 1e-15,
          reltol = 1e-15,
          p0 = nothing,
          verbose = false,
          threaded = true) where TT

  norm_prop,norm_property = normalize_property(model,prop,z,property)

  if norm_property !== property
    return _Pproperty(model,T,norm_prop,z,norm_property;rootsolver,phase,abstol,reltol,p0,verbose,threaded)
  end

  if length(model) == 1 && length(z) == 1
    res = Pproperty_pure(model,T,prop,z,property,rootsolver,phase,abstol,reltol,verbose,threaded,p0)
    return __Pproperty_check(res,verbose)
  end

  if p0 !== nothing
    res = __Pproperty(model,T,prop,z,property,rootsolver,phase,abstol,reltol,threaded,p0)
    return __Pproperty_check(res,verbose)
  end

  if is_liquid(phase)
    p00 = bubble_pressure(model,T,z)[1]
    res = __Pproperty(model,T,prop,z,property,rootsolver,phase,abstol,reltol,threaded,p00)
    return __Pproperty_check(res,verbose)
  end

  if is_vapour(phase)
    p00 = dew_pressure(model,T,z)[1]
    res = __Pproperty(model,T,prop,z,property,rootsolver,phase,abstol,reltol,threaded,p00)
    return __Pproperty_check(res,verbose)
  end

  n = sum(z)
  v0_edge,pure_sats = x0_edge_pressure(model,T,z)
  p0_bubble,p0_dew = v0_edge
  edge,crit,status = _edge_pressure(model,T,z,v0_edge)
  P_edge,v_l,v_v = edge

  if status == :supercritical
    #=
    TODO: what to do in this zone?
    we are smooth in the p-v curves,
    but there are still phase separation up until the mixture critical point. =#
    Tc,Pc,Vc = crit
    verbose && @info "mechanical critical pressure:        $Pc"
    verbose && @info "mechanical critical temperature:     $Tc"
    verbose && @info "mechanical critical molar volume     $Vc"
    res = __Pproperty(model,T,prop,z,property,rootsolver,:vapour,abstol,reltol,threaded,Pc)
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
      satpoint = "bubble"
      verbose && @info "molar volume at bubble point:        $Vsat"
    else
      dew_method_crit = dew_pressure_pproperty_method(model,Pc,Tc,z,pure_sats)
      Psat,_,Vsat,_ = dew_pressure(model,Tc,z,dew_method_crit)
      satpoint = "dew"
      verbose && @info "molar volume at dew point:           $Vsat"
    end

    verbose && @info "molar volume at pressure(property):  $Vx"

    βx = (Vx - Vsat)/(Vc - Vsat)
    0 <= βx <= 1 && verbose && @info "pseudo-critical pressure($property) in phase change region (between critical and $satpoint points)"
    0 <= βx <= 1 && return px,:eq
    verbose && @info "pressure(property) in the  critical pseudo-$(string(res[2])) branch, outside the phase change region"
    return res
  end

  if status == :failure
    verbose && @warn "failure to calculate edge point, trying to solve using Clapeyron.T_scale(model,z)"
    res = __Pproperty(model,T,prop,z,property,rootsolver,phase,abstol,reltol,threaded,p_scale(model,z))
    res[2] == :failure && return __Pproperty_check(res,verbose)
    return __Pproperty_check(res,verbose)
  end

  prop_l = spec_to_vt(model,v_l,T,z,property)
  prop_v = spec_to_vt(model,v_v,T,z,property)

  verbose && @info "property at liquid edge:     $prop_l"
  verbose && @info "property at vapour edge:     $prop_v"
  verbose && @info "pressure at edge point:      $P_edge"

  β = (prop - prop_l)/(prop_v - prop_l)

  #we are inside equilibria.
  if 0 <= β <= 1
    verbose && @info "property between the liquid and vapour edges, in the phase change region"
    return P_edge,:eq
  end

  new_phase = β > 1 ? :vapour : :liquid

  res = __Pproperty(model,T,prop,z,property,rootsolver,new_phase,abstol,reltol,threaded,P_edge)
  res[2] == :failure && return __Pproperty_check(res,verbose,P_edge)

  if β > 1 #check vapour branch
    dew_method = dew_pressure_pproperty_method(model,p0_dew,T,z,pure_sats)
    dew = dew_pressure(model,T,z,dew_method)
    p_dew,_,v_dew,_ = dew
    prob_dew = spec_to_vt(model,v_dew*n,T,z,property)

    verbose && @info "pressure at dew point:       $p_dew"
    verbose && @info "property at dew point:       $prob_dew"

    β_dew = (prop - prop_l)/(prob_dew - prop_l)

    0 < β_dew < 1 && verbose && @info "pseudo-vapour pressure($property) in phase change region (between edge and dew point)."
    0 < β_dew < 1 && return __Pproperty_check((res[1],:eq),verbose,P_edge)
    verbose && @info "vapour pressure($property) outside the phase change region"
    return __Pproperty_check(res,verbose)
  end

  if β < 0  #check liquid branch
    bubble_method = bubble_pressure_pproperty_method(model,p0_bubble,Tc,z,pure_sats)
    bubble = bubble_pressure(model,T,z,bubble_method)
    p_bubble,v_bubble,_,_ = bubble
    prob_bubble = spec_to_vt(model,v_bubble*n,T,z,property)

    verbose && @info "pressure at bubble point:    $p_bubble"
    verbose && @info "property at bubble point:    $prob_bubble"

    β_bubble = (prop - prop_l)/(prob_bubble - prop_l)

    0 < β_bubble < 1 && verbose && @info "pseudo-liquid pressure($property) in phase change region (between edge and bubble point)."
    0 < β_bubble < 1 && return __Pproperty_check((res[1],:eq),verbose,P_edge)
    verbose && @info "liquid pressure($property) outside the phase change region"
    return __Pproperty_check(res,verbose)
  end

  _0 = zero(Base.promote_eltype(model,T,prop,z))
  return __Pproperty_check((_0/_0,:failure),verbose)
end

function Pproperty_pure(model,T,x,z,property::F,rootsolver,phase,abstol,reltol,verbose,threaded,p0) where F
  TT = Base.promote_eltype(model,T,x,z)
  nan = zero(TT)/zero(TT)
  ∑z = sum(z)
  x1 = SA[1.0*one(∑z)]

  sat,crit,status = _extended_saturation_pressure(model,T)

  if status == :fail
    verbose && @error "PProperty calculation failed"
    return nan,:failure
  end

  if status == :supercritical
    verbose && @info "temperature is above critical temperature"
    Tc,Pc,Vc = crit
    return __Pproperty(model,T,x,z,property,rootsolver,phase,abstol,reltol,threaded,Pc)
  end

  ps,vl,vv = TT.(sat)

  xl = ∑z*spec_to_vt(model,vl,T,x1,property)
  xv = ∑z*spec_to_vt(model,vv,T,x1,property)
  βv = (x - xl)/(xv - xl)

  if !isfinite(βv)
    verbose && @error "PProperty calculation failed"
    return nan,:failure
  elseif βv < 0 || βv > 1
    phase0 = βv < 0 ? :liquid : :vapour
    is_liquid(phase0) && verbose && @info "pressure($property) > saturation pressure"
    is_vapour(phase0) && verbose && @info "pressure($property) < saturation pressure"
    return __Pproperty(model,T,x,z,property,rootsolver,phase0,abstol,reltol,threaded,ps)
  else
    verbose && @warn "$property value in phase change region. Will return pressure at saturation point"
    return ps,:eq
  end
end

function __Pproperty(model,T,prop,z,property::F,rootsolver,phase,abstol,reltol,threaded,p0) where F
  p_primal,phase = Pproperty_impl(primalval(model),primalval(T),primalval(prop),primalval(z),property,rootsolver,phase,abstol,reltol,threaded,primalval(p0))
  if has_a_res(model)
    p = Pproperty_ad(model,T,prop,z,property,p_primal,phase)
    return p,phase
  else
    p = p_primal
  end
  return p,phase
end

__Pproperty(model,T,prop,z,property::F,phase,p0) where F = __Pproperty(model,T,prop,z,property,Roots.Order0(),phase,1e-15,1e-15,true,p0)

function Pproperty_impl(model,T,prop,z,property::F,rootsolver,phase,abstol,reltol,threaded,p0) where F
  if is_unknown(phase)
    new_phase = identify_phase(model,p0,T,z)
    if is_unknown(new_phase) #something really bad happened
      _0 = zero(Base.promote_eltype(model,T,prop,z))
      nan = _0/_0
      return nan,:unknown
    end
    return __Pproperty(model,T,prop,z,property,rootsolver,new_phase,abstol,reltol,threaded,p0)
  end
  _1 = oneunit(typeof(prop))

  if property == volume
    px = pressure(model,prop,T,z)
    phasex = VT_identify_phase(model,prop,T,z)
    if is_unknown(phasex)
      return px,:failure
    else
      return px,phasex
    end
  end

  f(_lnp,_prop) = _1*property(model,exp(_lnp),T,z,phase = phase,threaded = threaded) - _prop
  prob = Roots.ZeroProblem(f,_1*log(p0))

  logp = Roots.solve(prob,rootsolver,p = prop,atol = abstol,rtol = reltol)
  if isnan(logp)
    return logp,:failure
  end
  return exp(logp),phase
end


function Pproperty_ad(model,T,prop,z,property::F,p_primal,phase) where F
  if has_dual(model) || has_dual(T) || has_dual(prop) || has_dual(z)
    #=
    we know that p_primal is the solution to
    property(model,p_primal,t,z,phase = phase,threaded = threaded) - prop = 0
    =#
    _property(_p) = property(model,_p,T,z,phase = phase)
    fprop,∂prop∂p = Solvers.f∂f(_property,p_primal)
    p = p_primal + (prop - fprop)/∂prop∂p
    return p
  else
    return p_primal
  end
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
