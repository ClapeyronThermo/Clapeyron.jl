function x0_Pproperty(model::EoSModel,T,z::AbstractVector,verbose = false)
  bubble = Clapeyron.bubble_pressure(model,T,z)
  dew = Clapeyron.dew_pressure(model,T,z)
  bubble_T = bubble[1]
  v_dew_vapour = dew[3]*sum(z)
  v_bubble_liquid = bubble[2]*sum(z)
  dew_T = dew[1]
  if isnan(bubble_T)
    verbose && @error "bubble_pressure calculation failed."
  end
  if isnan(dew_T)
    verbose && @error "dew_pressure calculation failed."
  end
  return (bubble_T,v_bubble_liquid),(dew_T,v_dew_vapour)
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
  if st == :failure
    @error "Pproperty calculation failed."
  end
  return p
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

  #handle volume variations
  if property == molar_density
    return _Pproperty(model,T,sum(z)/prop,z,volume;rootsolver,phase,abstol,reltol,p0,verbose,threaded)
  end

  if property == mass_density
    return _Pproperty(model,T,molecular_weight(model,z)/prop,z,volume;rootsolver,phase,abstol,reltol,p0,verbose,threaded)
  end

  if length(model) == 1 && length(z) == 1
    return Pproperty_pure(model,T,prop,z,property,rootsolver,phase,abstol,reltol,verbose,threaded,p0)
  end

  if p0 !== nothing
    return __Pproperty(model,T,prop,z,property,rootsolver,phase,abstol,reltol,threaded,p0)
  end

  if is_liquid(phase)
    p00 = bubble_pressure(model,T,z)[1]
    return __Pproperty(model,T,prop,z,property,rootsolver,phase,abstol,reltol,threaded,p00)
  end

  if is_vapour(phase)
    p00 = dew_pressure(model,T,z)[1]
    return __Pproperty(model,T,prop,z,property,rootsolver,phase,abstol,reltol,threaded,p00)
  end

  bubble_prop,dew_prop = x0_Pproperty(model,T,z,verbose)
  bubble_p,bubble_vol = bubble_prop
  dew_p,dew_vol = dew_prop

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
    verbose && @warn "non-finite bubble point, trying to solve using the dew point"
    return __Pproperty(model,T,prop,z,property,rootsolver,phase,abstol,reltol,threaded,dew_p)
  elseif !isnan(bubble_p) && isnan(dew_p)
    verbose && @warn "non-finite dew point, trying to solve using the bubble point"
    return __Pproperty(model,T,prop,z,property,rootsolver,phase,abstol,reltol,threaded,bubble_p)
  elseif isnan(bubble_p) && isnan(dew_p)
    verbose && @warn "non-finite dew and bubble points, trying to solve using Clapeyron.p_scale(model,z)"
    return __Pproperty(model,T,prop,z,property,rootsolver,phase,abstol,reltol,threaded,p_scale(model,z))
  end

  if property == volume
    prop_bubble = bubble_vol
    prop_dew = dew_vol
  else
    prop_bubble = property(model,bubble_p,T,z,phase=phase)
    prop_dew = property(model,dew_p,T,z,phase=phase)
  end

  F(P) = property(model,P,T,z,phase = phase)
  
  if verbose
    @info "input property:           $prop"
    @info "property at dew point:    $prop_dew"
    @info "property at bubble point: $prop_bubble"
    @info "pressure at dew point:    $dew_p"
    @info "pressure at bubble point: $bubble_p"
  end

  β = (prop - prop_dew)/(prop_bubble - prop_dew)
  if 0 <= β <= 1

      #special case with volume: the volumes in the saturation volumes don't have a definition of bulk:
      if property == volume
        verbose && @warn "volume in the phase change region, returning a linear interpolation of the bubble and dew pressures"
        return β*bubble_p + (1 - β)*dew_p,:eq
      end

      P_edge = FindEdge(F,dew_p,bubble_p) # dew_p < bubble_p --> condition for FindEdge
      if !isfinite(P_edge)
        verbose && @warn "failure to calculate edge point"
        verbose && @warn "$property in the phase change region, returning a linear interpolation of the bubble and dew pressures"
        return β*bubble_p + (1 - β)*dew_p,:eq
      end
      verbose && @info "pressure at edge point:   $P_edge"
      prop_edge1 = property(model,P_edge - 1e-10,T,z,phase = phase);
      prop_edge2 = property(model,P_edge + 1e-10,T,z,phase = phase);
      #=
      the order is the following:
      bubble -> edge1 -> edge2 -> dew
      or:
      dew -> edge2 -> edge1 -> bubble
      =#

      βedge = (prop - prop_edge1)/(prop_edge2 - prop_edge1)
   
      if 0 <= βedge <= 1
        verbose && @warn "In the phase change region"
        return P_edge,:eq
      elseif βedge < 0 #prop <= prop_edge2
        verbose && @info "pressure($property) ∈ (pressure(dew point),pressure(edge point))"
        px,_ = __Pproperty(model,T,prop,z,property,rootsolver,phase,abstol,reltol,threaded,dew_p)
        return px,:eq
        #abs(prop) > abs(prob_edge1)
      elseif βedge > 1
        verbose && @info "pressure($property) ∈ (pressure(edge point),pressure(bubble point))"
        px,_ = __Pproperty(model,T,prop,z,property,rootsolver,phase,abstol,reltol,threaded,bubble_p)
        return px,:eq
      end

    elseif β > 1
      verbose && @info "pressure($property) > pressure(bubble point)"
      return __Pproperty(model,T,prop,z,property,rootsolver,phase,abstol,reltol,threaded,bubble_p)
    elseif β < 0
      verbose && @info "pressure($property) < pressure(dew point)"
      return __Pproperty(model,T,prop,z,property,rootsolver,phase,abstol,reltol,threaded,dew_p)
    else
      _0 = Base.promote_eltype(model,T,prop,z)
      return _0/_0,:failure
    end
end

function Pproperty_pure(model,T,prop,z,property::F,rootsolver,phase,abstol,reltol,verbose,threaded,p0) where F
  ∑z = sum(z)
  if p0 !== nothing
    sol = __Pproperty(model,T,prop,z,property,rootsolver,phase,abstol,reltol,threaded,p0)
  end

  crit = crit_pure(model)
  Tc,Pc,Vc = crit

  if T >= Tc
    verbose && @info "temperature is above critical temperature"
    return __Pproperty(model,T,prop,z,property,rootsolver,phase,abstol,reltol,threaded,Pc)
  end

  Psat,vlsat,vvpat = saturation_pressure(model,T,crit = crit)

  if !is_unknown(phase)

    return __Pproperty(model,T,prop,z,property,rootsolver,phase,abstol,reltol,threaded,Psat)
  end

  if property == volume
    prop_v = vvsat
    prop_l = vlsat
  else
    prop_v = property(model,Psat,T,z,phase = :v)
    prop_l = property(model,Psat,T,z,phase = :l)
  end

  β = (prop - prop_l)/(prop_v - prop_l)
  if 0 <= β <= 1
    verbose && @warn "$property value in phase change region. Will return pressure at saturation point"
    return Psat,:eq
  elseif β < 0
    verbose && @info "pressure($property) > saturation pressure"
    return __Pproperty(model,T,prop,z,property,rootsolver,:liquid,abstol,reltol,threaded,Psat)

  elseif β > 1
    verbose && @info "pressure($property) < saturation pressure"
    return __Pproperty(model,T,prop,z,property,rootsolver,:vapour,abstol,reltol,threaded,Psat)
  else
    verbose && @error "PProperty calculation failed"
    _0 = Base.promote_eltype(model,T,prop,z)
    return _0/_0,:failure
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
  f(lnp,prop) = _1*property(model,exp(lnp),T,z,phase = phase,threaded = threaded) - prop
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
