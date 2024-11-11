"""
    `x0_Tproperty(model::EoSModel,p,z::AbstractVector)`
Peforms some initial checks to see if a possible solution exists in `Clapeyron.jl`.
"""
function x0_Tproperty(model::EoSModel,p,z::AbstractVector,verbose = false)
    @assert isapprox(sum(z),1,atol = 1e-8)
    bubble_prop = Clapeyron.bubble_temperature(model,p,z)
    dew_prop = Clapeyron.dew_temperature(model,p,z)
    bubble_temp = bubble_prop[1]
    dew_temp = dew_prop[1]
    if isnan(bubble_temp)
      verbose && @error "bubble_temp is NaN"
      return bubble_temp
    end
    if isnan(dew_temp)
      verbose && @error "dew_temp is NaN"
      return dew_temp
    end
    return bubble_temp,dew_temp
end

"""
    `FindEdge(f::Function,a,b)`
Finds approx singularity location in range `a`,`b` for function `f`. There should be only 1 singularity in [`a`,`b`].
"""
function FindEdge(f::Function,a,b)
  @assert b>= a
  if isapprox(a,b,atol=1e-10)
    return a
  end
    c = (a+b)/2
    f1,f2,f3 = f(a),f(c),f(b)
    ∇f1,∇f2 = (f2 - f1)/(c - a),(f3 - f2)/(b - a)
    if (∇f2 > ∇f1)
      FindEdge(f,c,b)
    else
      FindEdge(f,a,c)
    end
end

"""
    `Tproperty(model::EoSModel,p,prop,z::AbstractVector,property = enthalpy;rootsolver = Roots.Order0(),phase =:unknown,abstol = 1e-15,reltol = 1e-15, verbose = false)`

Given `p` and any other bulk property `prop` calculated via `property`, returns the required temperature `T` such that `property(model,p,T,z,phase) = prop`

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



  T,_ = _Tproperty(model,p,prop,z,property;rootsolver,phase,abstol,reltol,threaded,T0)
  return T
end

function _Tproperty(model::EoSModel,p,prop,z = SA[1.0],
                  property::TT = enthalpy;
                  rootsolver = Roots.Order0(),
                  phase =:unknown,
                  abstol = 1e-15,
                  reltol = 1e-15,
                  T0 = nothing,
                  verbose = false,
                  threaded = true) where TT


  if length(model) == 1 && length(z) == 1
    return Tproperty_pure(model,p,prop,property,rootsolver,phase,abstol,reltol,verbose,threaded,T0)
  end

  if T0 !== nothing
      return __Tproperty(model,p,prop,z,property,rootsolver,phase,abstol,reltol,threaded,T0)
  end

  bubble_temp,dew_temp = x0_Tproperty(model,p,z,verbose)

  #if any bubble/dew pressure is NaN, try solving for the non-NaN value
  #if both values are NaN, try solving using T_scale(model,z)
  if isnan(bubble_temp) && !isnan(dew_temp)
    verbose && @warn "non-finite bubble point, trying to solve using the dew point"
    return __Tproperty(model,p,prop,z,property,rootsolver,phase,abstol,reltol,threaded,dew_temp)
  elseif !isnan(bubble_temp) && isnan(dew_temp)
    verbose && @warn "non-finite dew point, trying to solve using the bubble point"
    return __Tproperty(model,p,prop,z,property,rootsolver,phase,abstol,reltol,threaded,bubble_temp)
  elseif isnan(bubble_temp) && isnan(dew_temp)
    verbose && @warn "non-finite dew and bubble points, trying to solve using Clapeyron.T_scale(model,z)"
    return __Tproperty(model,p,prop,z,property,rootsolver,phase,abstol,reltol,threaded,T_scale(model,z))
  end

  prop_bubble = property(model,p,bubble_temp,z,phase=phase)
  prop_dew = property(model,p,dew_temp,z,phase=phase)

  #case 1: Monotonically increasing
  if (prop_bubble < prop_dew)
    verbose && @info "$property at bubble < $property at dew point."
    if (prop < prop_bubble)
      verbose && @info "$property < $property at bubbble point"
      return __Tproperty(model,p,prop,z,property,rootsolver,phase,abstol,reltol,threaded,bubble_temp)
    end

    if prop_bubble <= prop <= prop_dew
      verbose && @info "$property at bubble <= $property <= $property at dew"

      T_edge = FindEdge(F,bubble_temp,dew_temp)
      prop_edge1 = property(model,p,T_edge - 1e-10,z,phase = phase)
      prop_edge2 = property(model,p,T_edge + 1e-10,z,phase = phase)

      if prop_edge1 <= prop <= prop_edge2
        verbose && @warn "In the phase change region"
        return T_edge,:eq 
      end

      if prop < prop_edge1
        verbose && @info "$property at bubble point < $property < $property at edge point"
        return __Tproperty(model,p,prop,z,property,rootsolver,phase,abstol,reltol,threaded,bubble_temp)
      end

      if prop > prop_edge2
        verbose && @info "$property at edge point < $property < $property at dew point"
        return __Tproperty(model,p,prop,z,property,rootsolver,phase,abstol,reltol,threaded,dew_temp)
      end
    end

    if prop > prop_dew
      verbose && @info "$property at dew point < $property"
      return __Tproperty(model,p,prop,z,property,rootsolver,phase,abstol,reltol,threaded,dew_temp)
    end
  end

  #case 2: Monotonically decreasing
  if prop_bubble > prop_dew
    verbose && @info "$property at bubble > $property at dew point."
    if prop > prop_bubble
      verbose && @info "$property > $property at bubbble point"
      return __Tproperty(model,p,prop,z,property,rootsolver,phase,abstol,reltol,threaded,bubble_temp)
    end

    if prop_dew <= prop <= prop_bubble
      verbose && @info "$property at bubble >= $property >= $property at dew"

      T_edge = FindEdge(F,bubble_temp,dew_temp)
      prop_edge1 = property(model,p,T_edge - 1e-10,z,phase = phase)
      prop_edge2 = property(model,p,T_edge + 1e-10,z,phase = phase)

      if prop_edge2 <= prop <= prop_edge1
        verbose && @warn "In the phase change region"
        return T_edge,:eq

      elseif prop >= prop_edge1
        verbose && @info "$property at bubble point > $property > $property at edge point"
        return __Tproperty(model,p,prop,z,property,rootsolver,phase,abstol,reltol,threaded,bubble_temp)
      end

      if prop <= prop_edge2
        verbose && @info "$property at edge point < $property < $property at dew point"
        return __Tproperty(model,p,prop,z,property,rootsolver,phase,abstol,reltol,threaded,dew_temp)
      end
    end

    if prop < prop_dew
      verbose && @info "$property at dew point > $property"
      return __Tproperty(model,p,prop,z,property,rootsolver,phase,abstol,reltol,threaded,dew_temp)
    end
  end

  _0 = Base.promote_eltype(model,p,prop,z)
  return _0/_0,:failure
end

function Tproperty_pure(model,p,prop,property::F,rootsolver,phase,abstol,reltol,verbose,threaded,T0) where F

  if T0 !== nothing
    sol = __Tproperty(model,p,prop,property,rootsolver,phase,abstol,reltol,threaded,T0)
  end

  crit = crit_pure(model)
  Tc,Pc,Vc = crit
  Tsat,vlsat,vvpat = Clapeyron.x0_saturation_temperature(model,p,crit)
  if isnan(Tsat)
    Tsat = critical_tsat_extrapolation(model,p,Tc,Pc,Vc)
  end

  Tmin = Tsat - 1
  Tmax = min(Tsat + 1,Tc)
  prop_edge1 = property(model,p,Tmin,phase = phase)
  prop_edge2 = property(model,p,Tmax,phase = phase)

  #case 1: property inside saturation dome
  if (prop_edge1 <= prop <= prop_edge2) || (prop_edge2 <= prop <= prop_edge1)
    verbose && @warn "$property value in phase change region. Will return temperature at saturation point"
    return Tsat,:eq
  end

  #case 2: Monotonically increasing property value i.e. prop_edge1 < prop_edge2
  if (prop_edge1 < prop_edge2)
    if (prop < prop_edge1)
      verbose && @info "$property < $property at saturation point"
      return __Tproperty(model,p,prop,property,rootsolver,phase,abstol,reltol,threaded,Tmin)
    else# (prop_edge2 < prop)
      verbose && @info "$property > $property at saturation point"
      return __Tproperty(model,p,prop,property,rootsolver,phase,abstol,reltol,threaded,Tmax)
    end
  else
  #case 3: Monotonically decreasing property value i.e. prop_edge1 > prop_edge2
    if (prop < prop_edge2)
      verbose && @info "$property < $property at saturation point"
      return __Tproperty(model,p,prop,property,rootsolver,phase,abstol,reltol,threaded,Tmax)
    else (prop_edge1 < prop)
      verbose && @info "$property > $property at saturation point"
      return __Tproperty(model,p,prop,property,rootsolver,phase,abstol,reltol,threaded,Tmin)
    end
  end
  _0 = Base.promote_eltype(model,p,prop,1.0)
  return _0/_0,:failure
end

function __Tproperty(model,p,prop,z,property::F,rootsolver,phase,abstol,reltol,threaded,T0) where F
  if is_unknown(phase)
    new_phase = identify_phase(model,p,T0,z)
    if is_unknown(new_phase) #something really bad happened
      _0 = zero(Base.promote_eltype(model,p,prop,z))
      nan = _0/_0
      return nan,:unknown
    end
    return __Tproperty(model,p,prop,z,property,rootsolver,new_phase,abstol,reltol,threaded,T0)
  end
  f(t,prop) = property(model,p,t,z,phase = phase,threaded = threaded) - prop
  prob = Roots.ZeroProblem(f,T0)
  sol = Roots.solve(prob,rootsolver,p = prop,atol = abstol,rtol = reltol)
  return sol,phase
end

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

export Tproperty