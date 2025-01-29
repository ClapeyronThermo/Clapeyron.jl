"""
    `x0_Tproperty(model::EoSModel,p,z::AbstractVector)`
Peforms some initial checks to see if a possible solution exists in `Clapeyron.jl`.
"""
function x0_Tproperty(model::EoSModel,p,z::AbstractVector,verbose = false)
    bubble = Clapeyron.bubble_temperature(model,p,z)
    dew = Clapeyron.dew_temperature(model,p,z)
    bubble_T = bubble[1]
    v_dew_vapour = dew[3]*sum(z)
    v_bubble_liquid = bubble[2]*sum(z)
    dew_T = dew[1]
    if isnan(bubble_T)
      verbose && @error "bubble_temperature calculation failed."
    end
    if isnan(dew_T)
      verbose && @error "dew_temperature calculation failed."
    end
    return (bubble_T,v_bubble_liquid),(dew_T,v_dew_vapour)
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
    if abs(∇f2) > abs(∇f1)
      FindEdge(f,c,b)
    else
      FindEdge(f,a,c)
    end
end

"""
    Tproperty(model::EoSModel,p,prop,z::AbstractVector,property = enthalpy;rootsolver = Roots.Order0(),phase =:unknown,abstol = 1e-15,reltol = 1e-15, verbose = false)

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

  T,st = _Tproperty(model,p,prop,z,property;rootsolver,phase,abstol,reltol,verbose,threaded,T0)
  if st == :failure
    @error "Tproperty calculation failed."
  end
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

  #handle volume variations
  if property == molar_density
    return _Tproperty(model,p,sum(z)/prop,z,volume;rootsolver,phase,abstol,reltol,verbose,threaded,T0)
  end

  if property == mass_density
    return _Tproperty(model,p,molecular_weight(model,z)/prop,z,volume;rootsolver,phase,abstol,reltol,verbose,threaded,T0)
  end

  if length(model) == 1 && length(z) == 1
    zz = SA[z[1]]
    return Tproperty_pure(model,p,prop,zz,property,rootsolver,phase,abstol,reltol,verbose,threaded,T0)
  end

  if T0 !== nothing
      return __Tproperty(model,p,prop,z,property,rootsolver,phase,abstol,reltol,threaded,T0)
  end

  if is_liquid(phase)
    T00 = bubble_temperature(model,p,z)[1]
    return __Tproperty(model,p,prop,z,property,rootsolver,phase,abstol,reltol,threaded,T00)
  end

  if is_vapour(phase)
    T00 = dew_temperature(model,p,z)[1]
    return __Tproperty(model,p,prop,z,property,rootsolver,phase,abstol,reltol,threaded,T00)
  end

  bubble_prop,dew_prop = x0_Tproperty(model,p,z,verbose)
  bubble_temp,bubble_vol = bubble_prop
  dew_temp,dew_vol = dew_prop

  #trivial
  if property === temperature
    T = prop*one(bubble_temp)
    β = (T - dew_temp)/(bubble_temp - dew_temp)
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
    return T,_new_phase
  end

  #if any bubble/dew temp is NaN, try solving for the non-NaN value
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

  if property == volume
    prop_bubble = bubble_vol
    prop_dew = dew_vol
  else
    prop_bubble = property(model,p,bubble_temp,z,phase=phase)
    prop_dew = property(model,p,dew_temp,z,phase=phase)
  end

  F(T) = property(model,p,T,z,phase = phase)

  if verbose
    @info "input property:              $prop"
    @info "property at dew point:       $prop_dew"
    @info "property at bubble point:    $prop_bubble"
    @info "temperature at dew point:    $dew_temp"
    @info "temperature at bubble point: $bubble_temp"
  end

  β = (prop - prop_dew)/(prop_bubble - prop_dew)
  if 0 <= β <= 1
      #special case with volume: the volumes in the saturation volumes don't have a definition of bulk:
      if property == volume
        verbose && @warn "volume in the phase change region, returning a linear interpolation of the bubble and dew temperatures"
        return β*bubble_temp + (1 - β)*dew_temp,:eq
      end

      T_edge = FindEdge(F,bubble_temp,dew_temp)
      if !isfinite(T_edge)
        verbose && @error "failure to calculate edge point"
        verbose && @warn "$property in the phase change region, returning a linear interpolation of the bubble and dew temperatures"
        return β*bubble_temp + (1 - β)*dew_temp,:eq
      end
      verbose && @info "temperature at edge point:   $T_edge"
      prop_edge1 = property(model,p,T_edge - 1e-10,z,phase = phase)
      prop_edge2 = property(model,p,T_edge + 1e-10,z,phase = phase)
      #=
      the order is the following:
      bubble -> edge1 -> edge2 -> dew
      or:
      dew -> edge2 -> edge1 -> bubble
      =#

      βedge = (prop - prop_edge1)/(prop_edge2 - prop_edge1)

      if 0 <= βedge <= 1
        verbose && @warn "In the phase change region"
        return T_edge,:eq
      elseif βedge > 1 #prop <= prop_edge2
        verbose && @info "temperature($property) ∈ (temperature(dew point),temperature(edge point))"
        return __Tproperty(model,p,prop,z,property,rootsolver,phase,abstol,reltol,threaded,dew_temp)
      elseif βedge < 0
        verbose && @info "temperature($property) ∈ (temperature(edge point),temperature(bubble point))"
        return __Tproperty(model,p,prop,z,property,rootsolver,phase,abstol,reltol,threaded,bubble_temp)
      end

    elseif β > 1
      verbose && @info "temperature($property) < temperature(bubble point)"
      __Tproperty(model,p,prop,z,property,rootsolver,phase,abstol,reltol,threaded,bubble_temp)
    elseif β < 0
      verbose && @info "temperature($property) > temperature(dew point)"
      __Tproperty(model,p,prop,z,property,rootsolver,phase,abstol,reltol,threaded,dew_temp)
    else
      verbose && @error "TProperty calculation failed"
      _0 = Base.promote_eltype(model,p,prop,z)
      return _0/_0,:failure
    end
end

function Tproperty_pure(model,p,prop,z,property::F,rootsolver,phase,abstol,reltol,verbose,threaded,T0) where F
  ∑z = sum(z)
  if T0 !== nothing
    sol = __Tproperty(model,p,prop,z,property,rootsolver,phase,abstol,reltol,threaded,T0)
  end

  crit = crit_pure(model)
  Tc,Pc,Vc = crit

  if p >= Pc
    verbose && @info "pressure is above critical pressure"
    return __Tproperty(model,p,prop,z,property,rootsolver,phase,abstol,reltol,threaded,Tc)
  end

  Tsat,vlsat,vvpat = saturation_temperature(model,p,crit = crit)

  if !is_unknown(phase)
    return __Tproperty(model,p,prop,z,property,rootsolver,phase,abstol,reltol,threaded,Tsat)
  end

  if property == volume
    prop_v = vvsat
    prop_l = vlsat
  else
    prop_v = property(model,p,Tsat,z,phase = :v)
    prop_l = property(model,p,Tsat,z,phase = :l)
  end

  β = (prop - prop_l)/(prop_v - prop_l)
  if 0 <= β <= 1
    verbose && @warn "$property value in phase change region. Will return temperature at saturation point"
    return Psat,:eq
  elseif β < 0
    verbose && @info "temperature($property) < saturation temperature"
    return __Tproperty(model,p,prop,z,property,rootsolver,:liquid,abstol,reltol,threaded,Tsat)

  elseif β > 1
    verbose && @info "temperature($property) > saturation temperature"
    return __Tproperty(model,p,prop,z,property,rootsolver,:vapour,abstol,reltol,threaded,Tsat)
  else
    verbose && @error "TProperty calculation failed"
    _0 = Base.promote_eltype(model,p,prop,z)
    return _0/_0,:failure
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
    return __Tproperty(model,p,prop,z,property,rootsolver,new_phase,abstol,reltol,threaded,T0)
  end
  f(t,prop) = property(model,p,t,z,phase = phase,threaded = threaded) - prop
  prob = Roots.ZeroProblem(f,T0)
  sol = Roots.solve(prob,rootsolver,p = prop,atol = abstol,rtol = reltol)
  if !isfinite(sol) || sol < 0
    return sol,:failure
  end
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
