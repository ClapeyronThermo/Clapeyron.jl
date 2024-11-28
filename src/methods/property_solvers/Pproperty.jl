function x0_Pproperty(model::EoSModel,T,z::AbstractVector,verbose = false)
  bubble = Clapeyron.bubble_pressure(model,T,z)
  dew = Clapeyron.dew_pressure(model,T,z)
  bubble_T = bubble[1]
  v_dew_vapour = dew[3]
  v_bubble_liquid = bubble[2]
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


  bubble_prop,dew_prop = x0_Pproperty(model,T,z,verbose)

  bubble_p,bubble_vol = bubble_prop
  dew_p,dew_vol = dew_prop

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

  #special case with volume: the volumes in the saturation volumes don't have a definition of bulk:
  if property == volume
    β = (prop - prop_dew)/(prop_bubble - prop_dew)
    if 0 <= β <= 1
      verbose && @warn "volume in the phase change region, returning a linear interpolation of the bubble and dew pressures"
      return β*bubble_p + (1 - β)*dew_p,:eq
    end
  end

  #case 1: Monotonically increasing
  if (prop_dew < prop_bubble)
    verbose && @info "$property at bubble > $property at dew point."
    if (prop < prop_dew)
      verbose && @info "$property < $property at dew point"
      return __Pproperty(model,T,prop,z,property,rootsolver,phase,abstol,reltol,threaded,dew_p)
    end

    if (prop >= prop_dew && prop <= prop_bubble)
      verbose && @info "$property at dew <= $property <= $property at bubble"
      P_edge = FindEdge(F,dew_p,bubble_p) # dew_p < bubble_pressure --> condition for FindEdge
      prop_edge1 = property(model,P_edge - 1e-10,T,z,phase = phase);
      prop_edge2 = property(model,P_edge + 1e-10,T,z,phase = phase);
      if (prop >= prop_edge1 && prop <= prop_edge2)
        verbose && @warn "In the phase change region"
        return P_edge,:eq
      end

      if (prop < prop_edge1)
        verbose && @info "$property at dew point < $property < $property at edge point"
        return __Pproperty(model,T,prop,z,property,rootsolver,phase,abstol,reltol,threaded,dew_p)
      end

      if (prop > prop_edge2)
        verbose && @info "$property at edge point < $property < $property at bubble point"
        return __Pproperty(model,T,prop,z,property,rootsolver,phase,abstol,reltol,threaded,bubble_p)
      end
    end

    if (prop > prop_bubble)
      verbose && @info "$property at bubble point < $property"
      __Pproperty(model,T,prop,z,property,rootsolver,phase,abstol,reltol,threaded,bubble_p)
    end
  end


  #case 2: Monotonically decreasing
  if (prop_bubble < prop_dew)
    verbose && @info "$property at bubble < $property at dew point."
    if (prop > prop_dew)
      verbose && @info "$property > $property at dew point"
      __Pproperty(model,T,prop,z,property,rootsolver,phase,abstol,reltol,threaded,dew_p)
    end

    if (prop<= prop_dew && prop >= prop_bubble)
      verbose && @info "$property at bubble <= $property <= $property at dew"
      P_edge = FindEdge(F,dew_p,bubble_p)

      prop_edge1 = property(model,P_edge - 1e-10,T,z,phase = phase)
      prop_edge2 = property(model,P_edge + 1e-10,T,z,phase = phase)
      @show prop_edge1,prop_edge2
      if (prop <= prop_edge1 && prop >= prop_edge2)
        verbose && @warn "In the phase change region"
        return P_edge,:eq

      elseif (prop >= prop_edge1)
        verbose && @info "$property at dew point > $property >= $property at edge point"
        return __Pproperty(model,T,prop,z,property,rootsolver,phase,abstol,reltol,threaded,dew_p)
      end

      if (prop <= prop_edge2)
        verbose && @info "$property at edge point >= $property > $property at bubble point"
        return __Pproperty(model,T,prop,z,property,rootsolver,phase,abstol,reltol,threaded,bubble_p)
      end
    end

    if (prop < prop_bubble)
      verbose && @info "$property at bubble point > $property"
      return __Pproperty(model,T,prop,z,property,rootsolver,phase,abstol,reltol,threaded,bubble_p)
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


