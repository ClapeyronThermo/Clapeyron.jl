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
    return (bubble_T,v_bubble_liquid),(dew_T,v_dew_vapour),(bubble,dew)
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


normalize_property(model,prop,z,property::F) where F = prop,property
normalize_property(model,prop,z,property::typeof(molar_density)) = sum(z)/prop,volume
normalize_property(model,prop,z,property::typeof(mass_density)) = molecular_weight(model,z)/prop,volume


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
  return T
end

function __Tproperty_check(res,verbose)
  T,st = res
  if verbose && st == :failure
    @error "TProperty calculation failed"
    return zero(T)/zero(T),st
  end
  return T,st
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


  norm_prop,norm_property = normalize_property(model,prop,z,property)

  if norm_property !== property
    res = _Tproperty(model,p,norm_prop,z,norm_property;rootsolver,phase,abstol,reltol,T0,verbose,threaded)
    return __Tproperty_check(res,verbose)
  end

  if length(model) == 1 && length(z) == 1
    zz = SA[z[1]]
    res = Tproperty_pure(model,p,prop,zz,property,rootsolver,phase,abstol,reltol,verbose,threaded,T0)
    return __Tproperty_check(res,verbose)
  end

  if T0 !== nothing
    res = __Tproperty(model,p,prop,z,property,rootsolver,phase,abstol,reltol,threaded,T0)
    return __Tproperty_check(res,verbose)
  end

  if is_liquid(phase)
    T00 = bubble_temperature(model,p,z)[1]
    res = __Tproperty(model,p,prop,z,property,rootsolver,phase,abstol,reltol,threaded,T00)
    return __Tproperty_check(res,verbose)
  end

  if is_vapour(phase)
    T00 = dew_temperature(model,p,z)[1]
    res = __Tproperty(model,p,prop,z,property,rootsolver,phase,abstol,reltol,threaded,T00)
    return __Tproperty_check(res,verbose)
  end

  bubble_prop,dew_prop,full_prop = x0_Tproperty(model,p,z,verbose)
  bubble_T,bubble_vol = bubble_prop
  dew_T,dew_vol = dew_prop

  #trivial
  if property === temperature
    T = prop*one(bubble_T)
    β = (T - dew_T)/(bubble_T - dew_T)
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
  if isnan(bubble_T) && !isnan(dew_T)
    verbose && @warn "non-finite bubble point, trying to solve using the dew point"
    _,dew_point = full_prop
    _,vl,vv,wdew = dew_point
    if property == volume
      prop_bubble = vl
      prop_dew = vv
    else
      prop_bubble = spec_to_vt(model,vl,T,wdew,spec)
      prop_dew = spec_to_vt(model,vv,T,z,spec)/sum(z)
    end
    β = (prop/sum(z) - prop_dew)/(prop_bubble - prop_dew)
    0 <= β <= 1 && return (dew_T,:eq)
    verbose && @info "pressure($property) < pressure(dew point)"
    return __Pproperty(model,p,prop,z,property,rootsolver,phase,abstol,reltol,threaded,dew_T)
  elseif !isnan(bubble_T) && isnan(dew_T)
    verbose && @warn "non-finite dew point, trying to solve using the bubble point"
        verbose && @warn "non-finite bubble point, trying to solve using the dew point"
        bubble_point,_ = full_prop
    _,vl,vv,wbubble = bubble_point
    if property == volume
      prop_bubble = vl
      prop_dew = vv
    else
      prop_bubble = spec_to_vt(model,vl,T,z,spec)/sum(z)
      prop_dew = spec_to_vt(model,vv,T,wbubble,spec)
    end

    β = (prop/sum(z) - prop_dew)/(prop_bubble - prop_dew)
    0 <= β <= 1 && (return bubble_T,:eq)
    verbose && @info "pressure($property) > pressure(bubble point)"
    return __Tproperty(model,p,prop,z,property,rootsolver,phase,abstol,reltol,threaded,bubble_T)
  elseif isnan(bubble_T) && isnan(dew_T)
    verbose && @warn "non-finite dew and bubble points, trying to solve using Clapeyron.p_scale(model,z)"
    return __Tproperty(model,p,prop,z,property,rootsolver,phase,abstol,reltol,threaded,p_scale(model,z))
  end

  if property == volume
    prop_bubble = bubble_vol
    prop_dew = dew_vol
  else
    prop_bubble = property(model,p,bubble_T,z,phase=phase)
    prop_dew = property(model,p,dew_T,z,phase=phase)
  end

  F(T) = property(model,p,T,z)

  if verbose
    @info "input property:              $prop"
    @info "property at dew point:       $prop_dew"
    @info "property at bubble point:    $prop_bubble"
    @info "temperature at dew point:    $dew_T"
    @info "temperature at bubble point: $bubble_T"
  end

  β = (prop - prop_dew)/(prop_bubble - prop_dew)

  if β > 1
    verbose && @info "temperature($property) < temperature(bubble point)"
    res = __Tproperty(model,p,prop,z,property,rootsolver,phase,abstol,reltol,threaded,bubble_T)
    return __Tproperty_check(res,verbose)
  elseif β < 0
    verbose && @info "temperature($property) > temperature(dew point)"
    res = __Tproperty(model,p,prop,z,property,rootsolver,phase,abstol,reltol,threaded,dew_T)
    return __Tproperty_check(res,verbose)
  end

  if 0 <= β <= 1
    T_edge = FindEdge(F,bubble_T,dew_T)
    if !isfinite(T_edge)
      verbose && @error "failure to calculate edge point"
      verbose && @warn "$property in the phase change region, returning a linear interpolation of the bubble and dew temperatures"
      return β*bubble_T + (1 - β)*dew_T,:eq
    end
    prop_edge_bubble = property(model,p,T_edge - 1e-10,z)
    prop_edge_dew = property(model,p,T_edge + 1e-10,z)
    verbose && @info "property at dew edge:        $prop_edge_dew"
    verbose && @info "property at bubble edge:     $prop_edge_bubble"
    verbose && @info "temperature at edge point:   $T_edge"

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
      return T_edge,:eq
    elseif βedge > 1
      verbose && @info "temperature($property) ∈ (temperature(dew point),temperature(edge point))"
      return βedge_dew*dew_T + (1 - βedge_dew)*T_edge,:eq
    elseif βedge < 0
      verbose && @info "temperature($property) ∈ (temperature(edge point),temperature(bubble point))"
      return βedge_bubble*bubble_T + (1 - βedge_bubble)*T_edge,:eq
    end
  end

  verbose && @error "TProperty calculation failed"
  _0 = zero(Base.promote_eltype(model,p,prop,z))
  return _0/_0,:failure
end

function Tproperty_pure(model,p,x,z,property::F,rootsolver,phase,abstol,reltol,verbose,threaded,T0) where F
    TT = Base.promote_eltype(model,p,x,z)
    nan = zero(TT)/zero(TT)
    ∑z = sum(z)
    x1 = SVector(1.0*one(∑z))

    sat,crit,status = _extended_saturation_temperature(model,p)

    if status == :fail
      verbose && @error "TProperty calculation failed"
      return nan,:failure
    end

    if status == :supercritical
      verbose && @info "pressure is above critical pressure"
      Tc,Pc,Vc = crit
      return __Tproperty(model,p,x,z,property,rootsolver,phase,abstol,reltol,threaded,Tc)
    end

    Ts,vl,vv = TT.(sat)

    xl = ∑z*spec_to_vt(model,vl,Ts,x1,property)
    xv = ∑z*spec_to_vt(model,vv,Ts,x1,property)
    βv = (x - xl)/(xv - xl)

    if !isfinite(βv)
      verbose && @error "TProperty calculation failed"
      return nan,:failure
    elseif βv < 0 || βv > 1
      phase0 = βv < 0 ? :liquid : :vapour
      is_liquid(phase0) && verbose && @info "temperature($property) < saturation temperature"
      is_vapour(phase0) && verbose && @info "temperature($property) > saturation temperature"
      return __Tproperty(model,p,x,z,property,rootsolver,phase0,abstol,reltol,threaded,Ts)
    else
      #verbose && @warn "$property value in phase change region. Will return temperature at saturation point"
      return Ts,:eq
    end
end

function __Tproperty(model,p,prop,z,property::F,rootsolver,phase,abstol,reltol,threaded,T0) where F
  T_primal,phase = Tproperty_impl(primalval(model),primalval(p),primalval(prop),primalval(z),property,rootsolver,phase,abstol,reltol,threaded,primalval(T0))
  if has_a_res(model)
    T = Tproperty_ad(model,p,prop,z,property,T_primal,phase)
    return T,phase
  else
    T = T_primal
  end
  return T,phase
end

__Tproperty(model,p,prop,z,property::F,phase,T0) where F = __Tproperty(model,p,prop,z,property,Roots.Order0(),phase,1e-15,1e-15,true,T0)

function Tproperty_impl(model,p,prop,z,property::F,rootsolver,phase,abstol,reltol,threaded,T0) where F
  if is_unknown(phase)
    new_phase = identify_phase(model,p,T0,z)
    if is_unknown(new_phase) #something really bad happened
      _0 = zero(Base.promote_eltype(model,p,prop,z))
      nan = _0/_0
      return nan,:unknown
    end
    return __Tproperty(model,p,prop,z,property,rootsolver,new_phase,abstol,reltol,threaded,T0)
  end
  _1 = oneunit(typeof(prop))
  f(t,prop) = _1*property(model,p,t,z,phase = phase,threaded = threaded) - prop
  prob = Roots.ZeroProblem(f,_1*T0)
  T = Roots.solve(prob,rootsolver,p = prop,atol = abstol,rtol = reltol)
  if !isfinite(T) || T < 0
    return T,:failure
  end
  return T,phase
end

function __Tproperty(model,p,prop,property::F,rootsolver,phase,abstol,reltol,threaded,T0) where F
  __Tproperty(model,p,prop,SA[1.0],property,rootsolver,phase,abstol,reltol,threaded,T0)
end

function Tproperty_ad(model,p,prop,z,property::F,T_primal,phase) where F
  if has_dual(model) || has_dual(p) || has_dual(prop) || has_dual(z)
    #=
    we know that T_primal is the solution to
    property(model,p,T_primal,z,phase = phase,threaded = threaded) - prop = 0
    =#
    _property(_T) = property(model,p,_T,z,phase = phase)
    fprop,∂prop∂T = Solvers.f∂f(_property,T_primal)
    T = T_primal + (prop - fprop)/∂prop∂T
    return T
  else
    return T_primal
  end
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
