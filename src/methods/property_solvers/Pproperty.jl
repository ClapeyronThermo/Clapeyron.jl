function x0_Pproperty(model::EoSModel,T,z::AbstractVector,verbose = false)
  bubble = Clapeyron.bubble_pressure(model,T,z)
  dew = Clapeyron.dew_pressure(model,T,z)
  if isnan(bubble[1])
    verbose && @error "bubble_pressure calculation failed."
  end
  if isnan(dew[1])
    verbose && @error "dew_pressure calculation failed."
  end
  return bubble,dew
end


"""
    Pproperty(model::EoSModel,T,prop,z = SA[1.0],property::TT = enthalpy;rootsolver = Roots.Order0(),phase =:unknown,abstol = 1e-15,reltol = 1e-15, verbose = false)

Given `T` and any other bulk property `prop` calculated via `property`, returns the required pressure `P` such that `property(model,p,T,z,phase) = prop`

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

  bubble,dew = x0_Pproperty(model,T,z,verbose)
  bubble_p,bubble_vl,bubble_vv,w_bubble = bubble
  dew_p,dew_vl,dew_vv,w_dew = dew

  #trivial
  if property === pressure
    p = prop*one(bubble_p)
    β = (p - bubble_p)/(dew_p - bubble_p)
    if 0 <= β <= 1
      verbose && @warn "In the phase change region"
      _new_phase = :eq
    elseif β > 1
      _new_phase = :vapour
    elseif β < 0
      _new_phase = :liquid
    else
      _new_phase = :failure
    end
    return p,_new_phase
  end

  #if any bubble/dew pressure is NaN, try solving for the non-NaN value
  #if both values are NaN, try solving using p_scale(model,z)
  if isnan(bubble_p) && !isnan(dew_p)

    #=
    this happens when there is a mixture with non-condensables (water + hydrogen)
    we suppose:
    - the properties of the liquid phase are not affected by changes in pressure (we solve for the gas phase)
    - the composition of the liquid is approximately equal to the composition of the liquid at the dew point.
    =#
  
    verbose && @warn "non-finite bubble point, trying to solve using the dew point"
    if property == volume
      prop_bubble = dew_vl
      prop_dew = dew_vv
    else
      prop_bubble = spec_to_vt(model,dew_vl,T,w_dew,property)
      prop_dew = spec_to_vt(model,dew_vv,T,z,property)/sum(z)
    end
    β = (prop/sum(z) - prop_bubble)/(prop_dew - prop_bubble)
    if 0 <= β <= 1
      #strategy: substract the liquid part, and solve for the gas fraction
      damp = Solvers.positive_linesearch(z/sum(z),w_dew,decay = 0.95) #make sure that the gas fraction is positive
      βx = damp*(1 - β)*sum(z)
      new_prop = prop - βx*prop_bubble
      new_z = z - w_dew * βx
      px,stx = __Pproperty(model,T,new_prop,new_z,property,rootsolver,:gas,abstol,reltol,threaded,dew_p)
      stx == :failure && (return dew_p,:eq)
      return px,:eq
    end
    verbose && @info "pressure($property) < pressure(dew point)"
    return __Pproperty(model,T,prop,z,property,rootsolver,phase,abstol,reltol,threaded,dew_p)
  elseif !isnan(bubble_p) && isnan(dew_p)

    #=
    this happens when there is a mixture with non-volatiles (water + dodecatriene)
    we suppose:
    - the properties of the liquid phase are not affected by changes in pressure (we solve for the gas phase)
    - the composition of the liquid is approximately equal to the composition of the liquid at the bubble point.
    =#

    verbose && @warn "non-finite dew point, trying to solve using the bubble point"
    if property == volume
      prop_bubble = bubble_vl
      prop_dew = bubble_vv
    else
      prop_bubble = spec_to_vt(model,bubble_vl,T,z,spec)/sum(z)
      prop_dew = spec_to_vt(model,bubble_vv,T,w_bubble,spec)
    end
    β = (prop/sum(z) - prop_dew)/(prop_bubble - prop_dew)
    if 0 <= β <= 1
      damp = Solvers.positive_linesearch(z/sum(z),w_bubble,decay = 0.95) #make sure that the gas fraction is positive
      βx = damp*(1 - β)*sum(z)
      new_prop = βx*prop_bubble
      new_z = w_bubble * βx
      px,stx = __Pproperty(model,T,new_prop,new_z,property,rootsolver,:gas,abstol,reltol,threaded,bubble_p)
      stx == :failure && (return bubble_p,:eq)
      return px,:eq
    end
    verbose && @info "pressure($property) > pressure(bubble point)"
    return __Pproperty(model,T,prop,z,property,rootsolver,phase,abstol,reltol,threaded,bubble_p)
  elseif isnan(bubble_p) && isnan(dew_p)
    verbose && @warn "non-finite dew and bubble points, trying to solve using Clapeyron.p_scale(model,z)"
    return __Pproperty(model,T,prop,z,property,rootsolver,phase,abstol,reltol,threaded,p_scale(model,z))
  end

  if property == volume
    prop_bubble = bubble_vl*sum(z)
    prop_dew = dew_vv*sum(z)
  else
    prop_bubble = property(model,bubble_p,T,z,phase=phase)
    prop_dew = property(model,dew_p,T,z,phase=phase)
  end

  F(P) = property(model,P,T,z)

  if verbose
    @info "input property:              $prop"
    @info "property at dew point:       $prop_dew"
    @info "property at bubble point:    $prop_bubble"
    @info "pressure at dew point:       $dew_p"
    @info "pressure at bubble point:    $bubble_p"
  end

  β = (prop - prop_dew)/(prop_bubble - prop_dew)

  if β > 1
    verbose && @info "pressure($property) > pressure(bubble point)"
    res = __Pproperty(model,T,prop,z,property,rootsolver,phase,abstol,reltol,threaded,bubble_p)
    return __Pproperty_check(res,verbose)
  elseif β < 0
    verbose && @info "pressure($property) < pressure(dew point)"
    res = __Pproperty(model,T,prop,z,property,rootsolver,phase,abstol,reltol,threaded,dew_p)
    return __Pproperty_check(res,verbose)
  elseif 0 <= β <= 1

    P_edge,prop_edge_dew,prop_edge_bubble = FindEdge(F,dew_p,bubble_p) # dew_p < bubble_p --> condition for FindEdge
    px0 = exp(β*log(bubble_p) + (1 - β)*log(dew_p))

    if !isfinite(P_edge)
      verbose && @warn "failure to calculate edge point"
      verbose && @warn "$property in the phase change region, returning a linear interpolation of the bubble and dew pressures"
      return px0,:eq
    end

    verbose && @info "property at dew edge:        $prop_edge_dew"
    verbose && @info "property at bubble edge:     $prop_edge_bubble"
    verbose && @info "pressure at edge point:      $P_edge"

    #=
    the order is the following:
    bubble -> edge_bubble -> edge_dew -> dew
    or:
    dew -> edge_dew -> edge_bubble -> bubble
    =#

    βedge = (prop - prop_edge_bubble)/(prop_edge_dew - prop_edge_bubble)
    βedge_bubble = (prop - prop_edge_bubble)/(prop_bubble - prop_edge_bubble)
    βedge_dew = (prop - prop_edge_dew)/(prop_dew - prop_edge_dew)

    if 0 <= βedge <= 1
      verbose && @warn "In the phase change region"
      return P_edge,:eq
    elseif βedge < 0 #prop <= prop_edge2
      verbose && @info "pressure($property) ∈ (pressure(bubble point),pressure(edge point))"
      P_edge_bubble = exp(βedge_bubble*log(bubble_p) + (1 - βedge_bubble)*log(P_edge))
      property == volume && return P_edge_bubble,:eq
      px,_ = __Pproperty(model,T,prop,z,property,rootsolver,phase,abstol,reltol,threaded,P_edge_bubble)
      return __Pproperty_check((px,:eq),verbose,P_edge_bubble)
    elseif βedge > 1
      verbose && @info "pressure($property) ∈ (pressure(edge point),pressure(dew point))"
      P_edge_dew = exp(βedge_dew*log(dew_p) + (1 - βedge_dew)*log(P_edge))
      property == volume && return P_edge_bubble,:eq
      px,_ = __Pproperty(model,T,prop,z,property,rootsolver,phase,abstol,reltol,threaded,P_edge_dew)
      return __Pproperty_check((px,:eq),verbose,P_edge_dew)
    end
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

export Pproperty
