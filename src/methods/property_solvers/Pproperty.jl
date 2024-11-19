function x0_Pproperty(model::EoSModel,T,z::AbstractVector,verbose = false)
    bubble_prop = Clapeyron.bubble_pressure(model,T,z)
    dew_prop = Clapeyron.dew_pressure(model,T,z)
    bubble_pressure = bubble_prop[1]
    dew_pressure = dew_prop[1]
    if isnan(bubble_pressure)
      verbose && @error "bubble_pressure is NaN"
      return bubble_pressure
    end
    if isnan(dew_pressure)
      verbose && @error "dew_pressure is NaN"
      return dew_pressure
    end
    return bubble_pressure,dew_pressure
end


"""
    ` Pproperty(model::EoSModel,T,prop,z = SA[1.0],property::TT = enthalpy;rootsolver = Roots.Order0(),phase =:unknown,abstol = 1e-15,reltol = 1e-15, verbose = false)`

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

  p,_ = _Pproperty(model,T,prop,z,property;rootsolver,phase,abstol,reltol,p0,verbose,threaded)
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

#   if length(model) == 1 && length(z) == 1
#     return Pproperty_pure(model::EoSModel,T,prop,property;rootsolver,phase,abstol,verbose)
#   end

  if p0 !== nothing
    return __Pproperty(model,T,prop,z,property,rootsolver,phase,abstol,reltol,threaded,p0)
  end

  bubble_pressure,dew_pressure = x0_Pproperty(model,T,z,verbose)
  #if any bubble/dew pressure is NaN, try solving for the non-NaN value
  #if both values are NaN, try solving using p_scale(model,z)
  if isnan(bubble_pressure) && !isnan(dew_pressure)
    verbose && @warn "non-finite bubble point, trying to solve using the dew point"
    return __Pproperty(model,T,prop,z,property,rootsolver,phase,abstol,reltol,threaded,dew_pressure)
  elseif !isnan(bubble_pressure) && isnan(dew_pressure)
    verbose && @warn "non-finite dew point, trying to solve using the bubble point"
    return __Pproperty(model,T,prop,z,property,rootsolver,phase,abstol,reltol,threaded,bubble_pressure)
  elseif isnan(bubble_pressure) && isnan(dew_pressure)
    verbose && @warn "non-finite dew and bubble points, trying to solve using Clapeyron.p_scale(model,z)"
    return __Pproperty(model,T,prop,z,property,rootsolver,phase,abstol,reltol,threaded,p_scale(model,z))
  end

  prop_bubble = property(model,bubble_pressure,T,z,phase=phase)
  prop_dew = property(model,dew_pressure,T,z,phase=phase)
  F(P) = property(model,P,T,z,phase = phase)
  #case 1: Monotonically increasing
  if (prop_dew < prop_bubble)
    verbose && @info "$property at bubble > $property at dew point."
    if (prop < prop_dew)
      verbose && @info "$property < $property at dew point"
      return __Pproperty(model,T,prop,z,property,rootsolver,phase,abstol,reltol,threaded,dew_pressure)
    end

    if (prop >= prop_dew && prop <= prop_bubble)
      verbose && @info "$property at dew <= $property <= $property at bubble"
      P_edge = FindEdge(F,dew_pressure,bubble_pressure) # dew_pressure < bubbble_pressure --> condition for FindEdge
      prop_edge1 = property(model,P_edge - 1e-10,T,z,phase = phase);
      prop_edge2 = property(model,P_edge + 1e-10,T,z,phase = phase);
      if (prop >= prop_edge1 && prop <= prop_edge2)
        verbose && @warn "In the phase change region"
        return P_edge,:eq
      end

      if (prop < prop_edge1)
        verbose && @info "$property at dew point < $property < $property at edge point"
        return __Pproperty(model,T,prop,z,property,rootsolver,phase,abstol,reltol,threaded,dew_pressure)
      end

      if (prop > prop_edge2)
        verbose && @info "$property at edge point < $property < $property at bubble point"
        return __Pproperty(model,T,prop,z,property,rootsolver,phase,abstol,reltol,threaded,bubble_pressure)
      end
    end

    if (prop > prop_bubble)
      verbose && @info "$property at bubble point < $property"
      __Pproperty(model,T,prop,z,property,rootsolver,phase,abstol,reltol,threaded,bubble_pressure)
    end
  end


  #case 2: Monotonically decreasing
  if (prop_bubble < prop_dew)
    verbose && @info "$property at bubble < $property at dew point."
    if (prop > prop_dew)
      verbose && @info "$property > $property at dew point"
      __Pproperty(model,T,prop,z,property,rootsolver,phase,abstol,reltol,threaded,dew_pressure)
    end

    if (prop<= prop_dew && prop >= prop_bubble)
      verbose && @info "$property at bubble <= $property <= $property at dew"
      P_edge = FindEdge(F,dew_pressure,bubble_pressure)

      prop_edge1 = property(model,P_edge - 1e-10,T,z,phase = phase)
      prop_edge2 = property(model,P_edge + 1e-10,T,z,phase = phase)
      if (prop <= prop_edge1 && prop >= prop_edge2)
        verbose && @warn "In the phase change region"
        return P_edge,:eq

      elseif (prop >= prop_edge1)
        verbose && @info "$property at dew point > $property >= $property at edge point"
        return __Pproperty(model,T,prop,z,property,rootsolver,phase,abstol,reltol,threaded,dew_pressure)
      end

      if (prop <= prop_edge2)
        verbose && @info "$property at edge point >= $property > $property at bubbble point"
        return __Pproperty(model,T,prop,z,property,rootsolver,phase,abstol,reltol,threaded,bubble_pressure)
      end
    end

    if (prop < prop_bubble)
      verbose && @info "$property at bubble point > $property"
      return __Pproperty(model,T,prop,z,property,rootsolver,phase,abstol,reltol,threaded,bubble_pressure)
    end
  end
  _0 = Base.promote_eltype(model,T,prop,z)
  return _0/_0,:failure
end

function Pproperty_pure(model,T,prop,z,property::F,rootsolver,phase,abstol,reltol,verbose,threaded,p0) where F
  ∑z = sum(z)
  if p0 !== nothing
    sol = __Pproperty(model,T,prop,z,property,rootsolver,phase,abstol,reltol,threaded,p0)
  end

  crit = crit_pure(model)
  Tc,Pc,Vc = crit
  Psat,vlsat,vvpat = saturation_pressure(model,T,crit = crit)
  if isnan(Tsat)
    psat = critical_psat_extrapolation(model,T,Tc,Pc,Vc)
  end

  Pmin = Psat*(1 - 1/400)
  Pmax = min(Psat*(1 + 1/400),Pc)
  prop_edge1 = property(model,Pmin,T,z,phase = phase)
  prop_edge2 = property(model,Pmax,T,z,phase = phase)

  #case 1: property inside saturation dome
  if T < Tc && (prop_edge1 <= prop <= prop_edge2) || (prop_edge2 <= prop <= prop_edge1)
    verbose && @warn "$property value in phase change region. Will return temperature at saturation point"
    return Psat,:eq
  end

  #case 2: Monotonically increasing property value i.e. prop_edge1 < prop_edge2
  if (prop_edge1 < prop_edge2)
    if (prop < prop_edge1)
      verbose && @info "$property < $property at saturation point"
      return __Pproperty(model,T,prop,z,property,rootsolver,phase,abstol,reltol,threaded,Pmin)
    else# (prop_edge2 < prop)
      verbose && @info "$property > $property at saturation point"
      return __Pproperty(model,T,prop,z,property,rootsolver,phase,abstol,reltol,threaded,Pmax)
    end
  else
  #case 3: Monotonically decreasing property value i.e. prop_edge1 > prop_edge2
    if (prop < prop_edge2)
      verbose && @info "$property < $property at saturation point"
      return __Pproperty(model,T,prop,z,property,rootsolver,phase,abstol,reltol,threaded,Pmax)
    else (prop_edge1 < prop)
      verbose && @info "$property > $property at saturation point"
      return __Pproperty(model,T,prop,z,property,rootsolver,phase,abstol,reltol,threaded,Pmin)
    end
  end
  _0 = Base.promote_eltype(model,T,prop,z)
  return _0/_0,:failure
end


function __Pproperty(model,T,prop,z,property::F,rootsolver,phase,abstol,reltol,threaded,p0) where F
  if is_unknown(phase)
    new_phase = identify_phase(model,p0,T,z)
    if is_unknown(new_phase) #something really bad happened
      _0 = zero(Base.promote_eltype(model,T,prop,z))
      nan = _0/_0
      return nan,:unknown
    end
    return __Pproperty(model,T,prop,z,property,rootsolver,new_phase,abstol,reltol,threaded,p0)
  end
  f(p,prop) = property(model,p,T,z,phase = phase,threaded = threaded) - prop
  prob = Roots.ZeroProblem(f,p0)
  sol = Roots.solve(prob,rootsolver,p = prop,atol = abstol,rtol = reltol)
  return sol,phase
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


