normalize_property(model,prop,z,property::F) where F = prop,property
normalize_property(model,prop,z,property::typeof(molar_density)) = sum(z)/prop,volume
normalize_property(model,prop,z,property::typeof(mass_density)) = molecular_weight(model,z)/prop,volume

function x0_edge_temperature(model,p,z,pure = split_pure_model(model))
  dPdTsat = extended_dpdT_temperature.(pure,p)
  T_bubble = antoine_bubble_solve(dPdTsat,p,z)
  T_dew = antoine_dew_solve(dPdTsat,p,z)
  return (T_bubble,T_dew),(pure,dPdTsat)
end

function μp_equality1_T2(model,p,z,x,Ts)
    lnv1,lnv2,T1,T2 = x
    n = sum(z)
    v1,v2 = exp(lnv1),exp(lnv2)
    RT1,RT2 = n*Rgas(model)*T1,n*Rgas(model)*T2
    f1(V) = a_res(model,V,T1,z)
    f2(V) = a_res(model,V,T2,z)
    A1,Av1 = Solvers.f∂f(f1,v1)
    A2,Av2 =Solvers.f∂f(f2,v2)
    p1,p2 = RT1*(-Av1 + 1/v1),RT2*(-Av2 + 1/v2)
    Δμᵣ = A1 - v1*Av1 - A2 + v2*Av2 + log(v2/v1)
    Fμ = Δμᵣ
    Fp1 = (p1 - p)/p
    Fp2 = (p2 - p)/p
    FT = (T1 - T2)/Ts
    return SVector(Fμ,Fp1,Fp2,FT)
end

mechanical_critical_point(model,z,x0) = crit_pure(model,x0,z)

function edge_temperature(model,p,z,v0 = nothing)
  edge,crit,status = _edge_temperature(model,p,z,v0)
  return edge
end

function _edge_temperature(model,p,z,v0 = nothing)
  if v0 == nothing
    vv0,_ = x0_edge_temperature(model,p,z)
  else
    vv0 = (v0[1],v0[2])
  end
  T1 = vv0[1]
  T2 = vv0[2]
  Tmin,Tmax = minmax(T1,T2)
  n = sum(z)
  v_Tmin = volume(model,p,Tmin,z,phase = :l)
  v_Tmax = volume(model,p,Tmax,z,phase = :v)
  Ts = 0.5*(T1 + T2)
  f(x) = μp_equality1_T2(model,p,z,x,Ts)
  V0 = SVector(promote(log(v_Tmin),log(v_Tmax),Tmin,Tmax))

  _0 = zero(V0[1])
  nan = _0/_0
  fail = (nan,nan,nan)

  _is_positive((v_Tmin,v_Tmax,Tmin,Tmax)) || return fail,fail,:failure

  sol = Solvers.nlsolve2(f,V0,Solvers.Newton2Var())
  v1 = exp(sol[1])
  v2 = exp(sol[2])
  T_eq = 0.5*(sol[3] + sol[4])
  edge = (T_eq,v1,v2)
  check_valid_sat_pure(model,p,v1,v2,T_eq,z) && (return edge,fail,:success)

  #fail when calculating edge temperature, this happens near the (mechanical) critical point
  Tr = T_eq/T_scale(model,z)
  vlog = log10(v1)
  crit = mechanical_critical_point(model,z,(Tr,vlog)) #mechanical critical point
  Tc,Pc,Vc = crit

  !isfinite(Pc) && return fail,fail,:failure
  Pc <= p && return fail,crit,:supercritical

  T_extrapolated = critical_tsat_extrapolation(model,p,Tc,Pc,Vc,z/sum(z))
  vlc,vvc = critical_vsat_extrapolation(model,T_extrapolated,Tc,Vc,z)
  V1 = SVector(promote(log(vlc),log(vvc),T_extrapolated,T_extrapolated))
  sol1 = Solvers.nlsolve2(f,V1,Solvers.Newton2Var())
  v3 = exp(sol1[1])
  v4 = exp(sol1[2])
  T_eq2 = 0.5*(sol1[3] + sol1[4])
  edge2 = (T_eq2,v3,v4)
  check_valid_sat_pure(model,p,v3,v4,T_eq2,z) && return edge2,crit,:success

  return fail,fail,:failure
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

function FindEdge(f::T,_a,_b,_fa,_fb) where T
  @assert _a <= _b

  a,b,fa,fb = promote(_a,_b,_fa,_fb)
  for i in 1:40
    isapprox(a,b,rtol=1e-10,atol = 1e-10) && return a,fa,fb
    c = 0.5*(a+b)
    fc = f(c)
    ∇fa,∇fc = (fc - fa)/(c - a),(fb - fc)/(b - a)
    if abs(∇fc) > abs(∇fa)
      a = c
      fa = fc
    else
      b = c
      fb = fc
    end
  end
  nan = zero(a)/zero(a)
  return nan,nan,nan
end

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

  edge,crit,status = _edge_temperature(model,p,z)
  T_edge,v_l,v_v = edge

  if status == :supercritical
    Tc,Pc,Vc = crit
    verbose && @info "mechanical critical pressure:        $Pc"
    verbose && @info "mechanical critical temperature:     $Tc"
    verbose && @info "temperature($property) over mechanical critical point"
    res = __Tproperty(model,p,prop,z,property,rootsolver,phase,abstol,reltol,threaded,Tc)
    res[2] == :failure && return __Tproperty_check(res,verbose)
    ψ_stable = diffusive_stability(model,p,res[1],z,phase = :vapour)
    !ψ_stable && verbose && @info "pseudo-critical temperature($property) in phase change region (diffusively unstable)"
    !ψ_stable && return __Tproperty_check((res[1],:eq),verbose)
    return __Tproperty_check(res,verbose)
  end

  if status == :failure
    verbose && @warn "failure to calculate edge point, trying to solve using Clapeyron.T_scale(model,z)"
    res = __Tproperty(model,p,prop,z,property,rootsolver,phase,abstol,reltol,threaded,T_scale(model,z))
    res[2] == :failure && return __Tproperty_check(res,verbose)
    ψ_stable = diffusive_stability(model,p,res[1],z,phase = :vapour)
    !ψ_stable && verbose && @info "pseudo-critical temperature($property) in phase change region (diffusively unstable)"
    !ψ_stable && return __Tproperty_check((res[1],:eq),verbose)
    return __Tproperty_check(res,verbose)
  end

  prop_l = spec_to_vt(model,v_l,T_edge,z,property)
  prop_v = spec_to_vt(model,v_v,T_edge,z,property)

  verbose && @info "property at liquid edge:     $prop_l"
  verbose && @info "property at vapour edge:     $prop_v"
  verbose && @info "temperature at edge point:   $T_edge"

  β = (prop - prop_l)/(prop_v - prop_l)

  #we are inside equilibria.
  if 0 <= β <= 1
    verbose && @info "property between the liquid and vapour edges, in the phase change region"
    return T_edge,:eq
  end

  #gas side, maybe eq, maybe not
  if β > 1
    res = __Tproperty(model,p,prop,z,property,rootsolver,:vapour,abstol,reltol,threaded,T_edge)
    res[2] == :failure && return __Tproperty_check(res,verbose,T_edge)
    ψ_stable = diffusive_stability(model,p,res[1],z,phase = :vapour)
    !ψ_stable && verbose && @info "pseudo-vapour temperature($property) in phase change region (diffusively unstable)"
    !ψ_stable && return __Tproperty_check((res[1],:eq),verbose,T_edge)

    verbose && @info "temperature($property) in vapour branch"

    dew = dew_temperature(model,p,z)
    T_dew,_,v_dew,_ = dew
    prob_dew = spec_to_vt(model,v_dew*sum(z),T_dew,z,property)
    
    verbose && @info "temperature at dew point:  $T_dew"
    verbose && @info "property at dew point:     $prob_dew"

    β_dew = (prop - prop_l)/(prob_dew - prop_l)
    0 < β_dew < 1 && verbose && @info "pseudo-liquid temperature($property) in phase change region (between edge and dew point)."
    0 < β_dew < 1 && return __Tproperty_check((res[1],:eq),verbose,T_edge)
    return __Tproperty_check(res,verbose)
  end

  if β < 0
    res = __Tproperty(model,p,prop,z,property,rootsolver,:liquid,abstol,reltol,threaded,T_edge)
    res[2] == :failure && return __Tproperty_check(res,verbose,T_edge)
    ψ_stable = diffusive_stability(model,p,res[1],z,phase = :liquid)

    !ψ_stable && verbose && @info "pseudo-liquid temperature($property) in phase change region (diffusively unstable)."
    !ψ_stable && return __Tproperty_check((res[1],:eq),verbose,T_edge)

    verbose && @info "temperature($property) in liquid branch"
    
    bubble = bubble_temperature(model,p,z)
    T_bubble,v_bubble,_,_ = bubble
    prob_bubble = spec_to_vt(model,v_bubble*sum(z),T_bubble,z,property)
    
    verbose && @info "temperature at bubble point:  $T_bubble"
    verbose && @info "property at bubble point:     $prob_bubble"

    β_bubble = (prop - prop_l)/(prob_bubble - prop_l)
    0 < β_bubble < 1 && verbose && @info "pseudo-liquid temperature($property) in phase change region (between edge and bubble point)."
    0 < β_bubble < 1 && return __Tproperty_check((res[1],:eq),verbose,T_edge)
    return __Tproperty_check(res,verbose)
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
