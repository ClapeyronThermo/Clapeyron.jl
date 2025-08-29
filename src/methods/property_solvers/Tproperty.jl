"""
    `x0_Tproperty(model::EoSModel,p,z::AbstractVector)`
Peforms some initial checks to see if a possible solution exists in `Clapeyron.jl`.
"""
function x0_Tproperty(model::EoSModel,p,z::AbstractVector,verbose = false)
    bubble = Clapeyron.bubble_temperature(model,p,z)
    dew = Clapeyron.dew_temperature(model,p,z)
    if isnan(bubble[1])
      verbose && @error "bubble_temperature calculation failed."
    end
    if isnan(dew[1])
      verbose && @error "dew_temperature calculation failed."
    end
    return bubble,dew
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

function FindEdge(f::T,a,b,fa,fb) where T
  @assert a <= b
  if isapprox(a,b,atol=1e-10)
    return a,fa,fb
  end
    c = (a+b)/2
    fc = f(c)
    ∇fa,∇fc = (fc - fa)/(c - a),(fb - fc)/(b - a)
    if abs(∇fc) > abs(∇fa)
      FindEdge(f,c,b,fc,fb)
    else
      FindEdge(f,a,c,fa,fc)
    end
end


normalize_property(model,prop,z,property::F) where F = prop,property
normalize_property(model,prop,z,property::typeof(molar_density)) = sum(z)/prop,volume
normalize_property(model,prop,z,property::typeof(mass_density)) = molecular_weight(model,z)/prop,volume


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
  T,st = _Tproperty(model,p,prop,z,property;rootsolver,phase,abstol,reltol,verbose,threaded,T0)
  return T
end

function __Tproperty_check(res,verbose,Tother = zero(res[1])/zero(res[1]))
  T,st = res
  if verbose && st == :failure
    @error "TProperty calculation failed"
    return Tother,st
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

  bubble,dew = x0_Tproperty(model,p,z,verbose)
  bubble_T,bubble_vl,bubble_vv,w_bubble = bubble
  dew_T,dew_vl,dew_vv,w_dew = dew

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
    if property == volume
      prop_bubble = dew_vl
      prop_dew = dew_vv
    else
      prop_bubble = spec_to_vt(model,dew_vl,dew_T,w_dew,property)
      prop_dew = spec_to_vt(model,dew_vv,dew_T,z,property)/sum(z)
    end
    β = (prop/sum(z) - prop_bubble)/(prop_dew - prop_bubble)
    if 0 <= β <= 1
      #strategy: substract the liquid part, and solve for the gas fraction
      #damp = Solvers.positive_linesearch(z/sum(z),w_dew,decay = 0.95) #make sure that the gas fraction is positive
      #βx = damp*(1 - β)*sum(z)
      #new_prop = prop - βx*prop_bubble
      #new_z = z - w_dew * βx
      #Tx,stx = __Tproperty(model,p,new_prop,new_z,property,rootsolver,:gas,abstol,reltol,threaded,dew_T)
      #stx == :failure && (return dew_T,:eq)
      return dew_T,:eq
    end
    verbose && @info "pressure($property) < pressure(dew point)"
    return __Pproperty(model,p,prop,z,property,rootsolver,phase,abstol,reltol,threaded,dew_T)
  elseif !isnan(bubble_T) && isnan(dew_T)
    verbose && @warn "non-finite dew point, trying to solve using the bubble point"
    if property == volume
      prop_bubble = bubble_vl
      prop_dew = bubble_vv
    else
      prop_bubble = spec_to_vt(model,bubble_vl,bubble_T,z,property)/sum(z)
      prop_dew = spec_to_vt(model,bubble_vv,bubble_T,w_bubble,property)
    end

    β = (prop/sum(z) - prop_dew)/(prop_bubble - prop_dew)
    0 <= β <= 1 && (return bubble_T,:eq)
    verbose && @info "pressure($property) > pressure(bubble point)"
    return __Tproperty(model,p,prop,z,property,rootsolver,phase,abstol,reltol,threaded,bubble_T)
  elseif isnan(bubble_T) && isnan(dew_T)
    verbose && @warn "non-finite dew and bubble points, trying to solve using Clapeyron.T_scale(model,z)"
    return __Tproperty(model,p,prop,z,property,rootsolver,phase,abstol,reltol,threaded,T_scale(model,z))
  end

  if property == volume
    prop_bubble = bubble_vl*sum(z)
    prop_dew = dew_vv*sum(z)
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
  elseif 0 <= β <= 1

    if property == volume
      verbose && @info "$property in the phase change region, returning a linear interpolation of the bubble and dew temperatures"
      return β*bubble_T + (1 - β)*dew_T,:eq
    end

    T_edge,prop_edge_bubble,prop_edge_dew = FindEdge(F,bubble_T,dew_T)

    if !isfinite(T_edge)
      verbose && @error "failure to calculate edge point"
      verbose && @warn "$property in the phase change region, returning a linear interpolation of the bubble and dew temperatures"
      return β*bubble_T + (1 - β)*dew_T,:eq
    end
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
      T_edge_dew = βedge_dew*dew_T + (1 - βedge_dew)*T_edge
      #TODO: we could skip this calculation when Tproperty is used as initial point.
      Tx,_ = __Tproperty(model,p,prop,z,property,rootsolver,phase,abstol,reltol,threaded,T_edge_dew)
      return __Tproperty_check((Tx,:eq),verbose,T_edge)
    elseif βedge < 0
      verbose && @info "temperature($property) ∈ (temperature(edge point),temperature(bubble point))"
      T_edge_bubble = βedge_bubble*bubble_T + (1 - βedge_bubble)*T_edge
      #TODO: we could skip this calculation when Tproperty is used as initial point.
      Tx,_ = __Tproperty(model,p,prop,z,property,rootsolver,phase,abstol,reltol,threaded,T_edge_bubble)
      return __Tproperty_check((Tx,:eq),verbose,T_edge)
    end
  end

  _0 = zero(Base.promote_eltype(model,p,prop,z))
  return __Tproperty_check((_0/_0,:failure),verbose)
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

  #return Tproperty_solver(model,p,prop,z,property,phase,abstol,reltol,T0)
end
#=
function Tproperty_solver(model,p,prop,z,property,phase = :unknown,abstol = 1e-15,reltol = 1e-15,T0 = NaN)
  XX = Base.promote_eltype(model,p,prop,z)
  Ta = XX(T0)
  h = cbrt(eps(one(Ta)))
  δT = h * oneunit(Ta) + abs(Ta) * h^2

  Tb = T0 + δT
  nan = (0Ta)/(0Ta)
  Tmin,Tmax = nan,nan
  vmin,vmax = nan,nan
  va::XX = volume(model,p,Ta,z,phase = phase)
  vb::XX = volume(model,p,Tb,z,phase = phase)
  fa::XX = spec_to_vt(model,va,Ta,z,property) - prop
  fb::XX = spec_to_vt(model,vb,Tb,z,property) - prop
  abs(fa) <= max(abstol, abs(Ta) * reltol) && return Ta,phase
  abs(fb) <= max(abstol, abs(Tb) * reltol) && return Tb,phase
  fa == fb && return nan,:failure

  #step 1: secant
  success = false
  for _ in 1:100
    Tm::XX = Tb - (Tb - Ta) * fb / (fb - fa)
    vm::XX = volume(model,p,Tm,z,phase = phase)
    fm::XX = spec_to_vt(model,vm,Tm,z,property) - prop
    iszero(fm) && return Tm,phase
    isnan(fm) && return nan,:failure
    abs(fm) <= max(abstol, abs(Tm) * reltol) && return Tm,phase
    if fm == fb
      return nan,phase
    end
    Tmin,Tmax = minmax(Ta,Tb)
    vmin,vmax = minmax(va,vb)
    if Tmin <= Tm <= Tmax
      success = true
      break
    end
    Ta, Tb, fa, fb, va, vb = Tb, Tm, fb, fm, vb, vm
  end
  success || (return nan,:failure)
  #step 2: newton
  f_newton(vt) = Tproperty_obj(vt,model,p,prop,z,property)
  fj(xx) = Solvers.FJ_ad(f_newton,xx)
  x = SVector(0.5*(vmin + vmax),0.5*(Tmin + Tmax))
  for _ in 1:20
    Fx,Jx = fj(x)
    d = Jx \ -Fx
    y = x + d
    y1,y2 = y
    y1 < vmin && (y1 = 0.5*(x[1] + vmin))
    y1 > vmax && (y1 = 0.5*(x[1] + vmax))
    y2 < Tmin && (y2 = 0.5*(x[2] + Tmin))
    y2 > Tmax && (y2 = 0.5*(x[2] + Tmax))
    x = SVector(y1,y2)
    ρF = norm(Fx, Inf)
    ρs = norm(d, Inf)
    ρx = norm(x, Inf)
    #@show ρF, ρs
    if ρs <= max(abstol, ρx*reltol) || ρF <= max(abstol, ρx * reltol)
        return x[2],phase
    end

    if !all(isfinite,x)
        return nan,:failure
    end
  end
  return nan,:failure
end

function Tproperty_obj(vt,model,p,x,z,spec)
  v,T = vt
  px = pressure(model,v,T,z)
  propx = spec_to_vt(model,v,T,z,spec)
  F1 = (p - px)/p
  F2 = (propx - x)/x
  return SVector(F1,F2)
end
=#
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
