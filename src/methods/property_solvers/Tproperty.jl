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
    `FindEdge(f::Function,a::Float64,b::Float64)`
Finds approx singularity location in range `a`,`b` for function `f`. There should be only 1 singularity in [`a`,`b`].
"""
function FindEdge(f::Function,a,b)
  @assert b>= a
  if isapprox(a,b,atol=1e-10)
    return a;
  end
    c = (a+b)/2;
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
function Tproperty(model::EoSModel,p,prop,z = SA[1.0],property::TT = enthalpy;rootsolver = Roots.Order0(),phase =:unknown,abstol = 1e-15,reltol = 1e-15, verbose = false) where TT
  if length(model) == 1 && length(z) == 1
    return Tproperty_pure(model::EoSModel,p,prop,property;rootsolver,phase,abstol,verbose)
  end
  bubble_temp,dew_temp = x0_Tproperty(model,p,z,verbose)
  if isnan(bubble_temp) || isnan(dew_temp)
    #TODO handle this case
    return bubble_temp
  end
  prop_bubble = property(model,p,bubble_temp,z,phase=phase)
  prop_dew = property(model,p,dew_temp,z,phase=phase)
  f(t,prop) = property(model,p,t,z,phase = phase) - prop
  F(T) = property(model,p,T,z,phase = phase)
  #case 1: Monotonically increasing
  if (prop_bubble < prop_dew)
    verbose && @info "$property at bubble < $property at dew point."
    if (prop < prop_bubble)
      verbose && @info "$property < $property at bubbble point"
      prob = Roots.ZeroProblem(f,bubble_temp)
      sol = Roots.solve(prob,rootsolver,prop,atol = abstol,rtol = reltol)
      return sol
    end

    if (prop >= prop_bubble && prop <= prop_dew)
      verbose && @info "$property at bubble <= $property <= $property at dew"
      T_edge = FindEdge(F,bubble_temp,dew_temp)
      prop_edge1 = property(model,p,T_edge - 1e-10,z,phase = phase);
      prop_edge2 = property(model,p,T_edge + 1e-10,z,phase = phase);
      if (prop >= prop_edge1 && prop <= prop_edge2)
        verbose && @warn "In the phase change region"
        return T_edge
      end

      if (prop < prop_edge1)
        verbose && @info "$property at bubble point < $property < $property at edge point"
        prob = Roots.ZeroProblem(f,bubble_temp)
        sol = Roots.solve(prob,rootsolver,prop,atol = abstol,rtol = reltol)
        return sol
      end

      if (prop>prop_edge2)
        verbose && @info "$property at edge point < $property < $property at dew point"
        prob = Roots.ZeroProblem(f,dew_temp)
        sol = Roots.solve(prob,rootsolver,prop,atol = abstol,rtol = reltol)
        return sol
      end
    end

    if (prop > prop_dew)
      verbose && @info "$property at dew point < $property"
      prob = Roots.ZeroProblem(f,dew_temp)
      sol = Roots.solve(prob,rootsolver,prop,atol = abstol,rtol = reltol)
      return sol
    end
  end


  #case 2: Monotonically decreasing
  if (prop_bubble > prop_dew)
    verbose && @info "$property at bubble > $property at dew point."
    if (prop > prop_bubble)
      verbose && @info "$property > $property at bubbble point"
      prob = Roots.ZeroProblem(f,bubble_temp)
      sol = Roots.solve(prob,rootsolver,prop,atol = abstol,rtol = reltol)
      return sol
    end

    if (prop<= prop_bubble && prop >= prop_dew)
      verbose && @info "$property at bubble >= $property >= $property at dew"
      T_edge = FindEdge(F,bubble_temp,dew_temp)

      prop_edge1 = property(model,p,T_edge - 1e-10,z,phase = phase)
      prop_edge2 = property(model,p,T_edge + 1e-10,z,phase = phase)
      if (prop <= prop_edge1 && prop >= prop_edge2)
        verbose && @warn "In the phase change region"
        return T_edge

      elseif (prop >= prop_edge1)
        verbose && @info "$property at bubble point > $property > $property at edge point"
        prob = Roots.ZeroProblem(f,bubble_temp)
        sol = Roots.solve(prob,rootsolver,prop,atol = abstol,rtol = reltol)
        return sol
      end

      if (prop <= prop_edge2)
        verbose && @info "$property at edge point < $property < $property at dew point"
        prob = Roots.ZeroProblem(f,dew_temp)
        sol = Roots.solve(prob,rootsolver,prop,atol = abstol,rtol = reltol)
        return sol
      end
    end

    if (prop < prop_dew)
      verbose && @info "$property at dew point > $property"
      prob = Roots.ZeroProblem(f,dew_temp)
      sol = Roots.solve(prob,rootsolver,prop,atol = abstol,rtol = reltol)
      return sol
    end
  end
end



"""
`Tproperty(model::EoSModel,p,prop,property::Function;rootsolver = FastShortcutNonlinearPolyalg(),phase =:unknown)` for single components `model`

Requires `p` and any other bulk property `prop` to compute the necessary temperature.
"""
function Tproperty_pure(model::EoSModel,p,prop,property::Function;rootsolver = Roots.Order0(),phase =:unknown,abstol = 1e-15,reltol = 1e-15,verbose=false)
  crit = crit_pure(model)
  Tc,Pc,Vc = crit
  Tsat,vlsat,vvpat = Clapeyron.x0_saturation_temperature(model,p,crit)
  if isnan(Tsat)
    if p < Pc
      verbose && @error "Saturation temperature not found"
      return Tsat
    else
      #use saturation extrapolation. we extrapolate the saturation curve based on the slope at the critical point and calculate a "pseudo" temperature.
      _p(_T) = pressure(pure,Vc,_T)
      dpdT = Solvers.derivative(_p,Tc)
      dlnPdTinvsat = -dpdT*Tc*Tc/Pc
      Δlnp = log(p/pii)
      #dT = clamp(dTdp*Δp,-0.5*T,0.5*T)
      Tinv0 = 1/T
      Tinv = Tinv0 + dTinvdlnp*Δlnp
      Tsat = 1/Tinv
    end
  end
  Tmin = Tsat - 1
  Tmax = min(Tsat + 1,Tc)
  prop_edge1 = property(model,p,Tmin,phase = phase)
  prop_edge2 = property(model,p,Tmax,phase = phase)
  f(t,prop) = property(model,p,t,phase = phase) - prop;
  #case 1: Monotonically increasing property value i.e. prop_edge1 < prop_edge2
  if (prop_edge1 < prop_edge2)
    if (prop < prop_edge1)
      verbose && @info "$property < $property at saturation point"
      prob = Roots.ZeroProblem(f,Tmin)
      sol = Roots.solve(prob,rootsolver,prop,atol = abstol,rtol = reltol)
    elseif (prop_edge1 <= prop && prop_edge2 >= prop)
      verbose && @warn "$property value in phase change region. Will return temperature at saturation point"
      sol = Tsat
    elseif (prop_edge2< prop)
      verbose && @info "$property > $property at saturation point"
      prob = Roots.ZeroProblem(f,Tmax)
      sol = Roots.solve(prob,rootsolver,prop,atol = abstol,rtol = reltol)
      return sol
    end
  end
  #case 2: Monotonically decreasing property value i.e. prop_edge1 > prop_edge2
  if (prop_edge2 < prop_edge1)
    if (prop < prop_edge2)
      verbose && @info "$property < $property at saturation point"
      prob = Roots.ZeroProblem(f,Tmax)
      sol = Roots.solve(prob,rootsolver,prop,atol = abstol,rtol = reltol)
    elseif (prop_edge1 >= prop && prop_edge2 <= prop)
      verbose && @warn "$property value in phase change region. Will return temperature at saturation point"
      sol = Tsat
    elseif (prop_edge1 < prop)
      verbose && @info "$property > $property at saturation point"
      prob = Roots.ZeroProblem(f,Tmin)
      sol = Roots.solve(prob,rootsolver,prop,atol = abstol,rtol = reltol)
    end
  end
  return sol
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