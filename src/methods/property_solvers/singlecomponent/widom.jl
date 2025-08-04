#=
widom line:
CpdT|p = 0
dCₚdT|p = dCₚdT|V - dCₚdV|T * dpdT|V /dpdV|T

CIIC line:
dCₚdP|T = 0
dCₚdP|T = dCₚdV|T * dPdV|T
=
=#
"""
    pw, vw = widom_pressure(model, T; p0 = nothing, v0 = nothing)

Calculate the Widom pressure and corresponding volume at a given temperature for a thermodynamic model.

The Widom point represents the temperature `T` where maxima of the isobaric heat capacity occurs along an isobar in the supercritical region.
In particular, `widom_pressure` calculates the widom point at a specified temperature.

# Arguments
- `model`: Thermodynamic model used for calculations.
- `T`: Temperature at which to calculate the Widom pressure `[K]`.
- `p0`: Optional initial pressure guess `[Pa]`. If not provided, a default value will be used.
- `v0`: Optional initial volume guess `[m³]`. If not provided, a default value will be used.

# Returns
- `pw`: Widom pressure at the specified temperature `[Pa]`.
- `vw`: Volume at the Widom pressure and specified temperature `[m³]`.

# Example
```julia
model = PCSAFT(["methane"])
T = 200.0  # 200 K
pw, vw = widom_pressure(model, T)
```
"""
function widom_pressure(model,T;p0 = nothing,v0 = nothing)
    if !isnothing(p0) && isnothing(v0)
        _v0 = volume(model,p0,T,phase = :l)
    elseif isnothing(p0) && !isnothing(v0)
        _v0 = v0*one(Base.promote_eltype(model,T))
    elseif isnothing(p0) && isnothing(v0)
        crit = crit_pure(model)
        p0 = critical_psat_extrapolation(model,T,crit)
        return widom_pressure(model,T;p0 = p0)
    else
        throw(ArgumentError("cannot specify `v0` and `p0` simultaneously"))
    end

    function f_wp(lnv)
        v = exp(lnv)
        Cₚ(∂v,∂T) = VT_isobaric_heat_capacity(model,∂v,∂T,SA[1.0])
        p(∂v,∂T) = pressure(model,∂v,∂T,SA[1.0])
        ∂Cₚ∂V,∂Cₚ∂T = Solvers.gradient2(Cₚ,v,T)
        ∂p∂V,∂p∂T = Solvers.gradient2(p,v,T)
        return ∂Cₚ∂T - ∂Cₚ∂V * ∂p∂T / ∂p∂V
    end

    log_v_widom = Solvers.nlsolve(f_wp,log(_v0),Roots.Order0())
    v_widom = exp(log_v_widom)
    p_widom = pressure(model,v_widom,T,SA[1.0])
    return p_widom,v_widom
end

"""
    Tw, vw = widom_temperature(model, p; T0 = nothing, v0 = nothing)

Calculate the Widom temperature and corresponding volume at a given pressure for a thermodynamic model.

The Widom point represents the temperature `T` where maxima of the isobaric heat capacity occurs along an isobar in the supercritical region.
In particular, `widom_temperature` calculates the widom point at a specified temperature.

# Arguments
- `model`: Thermodynamic model used for calculations.
- `p`: Pressure at which to calculate the Widom temperature `[Pa]`.
- `T0`: Optional initial temperature guess `[K]`. If not provided, a default value will be used.
- `v0`: Optional initial volume guess `[m³]`. If not provided, a default value will be used.

# Returns
- `Tw`: Widom temperature at the specified pressure `[K]`.
- `vw`: Volume at the Widom temperature and specified pressure `[m³]`.

# Example
```julia
model = PCSAFT(["methane"])
p = 150e5  #150 bar, over critical pressure
Tw, vw = widom_temperature(model, p)
```
"""
function widom_temperature(model,p;T0 = nothing,v0 = nothing)
    if !isnothing(T0) && isnothing(v0)
        _v0 = volume(model,p,T0,phase = :l)
        _T0 = T0*one(Base.promote_eltype(model,p))
    elseif isnothing(T0) && !isnothing(v0)
        _v0 = v0*one(Base.promote_eltype(model,p))
        crit = crit_pure(model)
        _T0 = critical_tsat_extrapolation(model,p,crit)
    elseif isnothing(T0) && isnothing(v0)
        crit = crit_pure(model)
        T0 = critical_tsat_extrapolation(model,p,crit)
        return widom_temperature(model,p;T0 = T0)
    else
        _1 = Base.promote_eltype(model,p,T0,v0)
        _v0 = v0*_1
        _T0 = T0*_1
    end

    function f_wp(x::AbstractVector)
        v,T = exp(x[1]),x[2]
        Cₚ(∂v,∂T) = VT_isobaric_heat_capacity(model,∂v,∂T,SA[1.0])
        fp(∂v,∂T) = pressure(model,∂v,∂T,SA[1.0])
        ∂Cₚ∂V,∂Cₚ∂T = Solvers.gradient2(Cₚ,v,T)
        px,dpx = Solvers.fgradf2(fp,v,T)
        ∂p∂V,∂p∂T = dpx
        widom_condition = (∂Cₚ∂T - ∂Cₚ∂V * ∂p∂T / ∂p∂V)#*T/Rgas(model)
        pressure_eq = (px - p)/p
        return SVector(widom_condition,pressure_eq)
    end

    #we do a refinement of the initial value first. it improves convergence near the critical point
    function f_wp(T)
        v = volume(model,p,T,phase = :l,vol0 = _v0)
        return f_wp(SVector(log(v),T))[1]
    end
    prob = Roots.ZeroProblem(f_wp,_T0)
    _T00 = Roots.solve(prob,Roots.Order2(),atol = 1e-1)
    _v00 = volume(model,p,_T00,phase = :l, vol0 = _v0)
    res = Solvers.nlsolve2(f_wp,SVector(log(_v00),_T00),Solvers.Newton2Var())
    log_v_widom,T_widom = res
    return T_widom,exp(log_v_widom)
end

"""
    p_ciic, v_ciic = ciic_pressure(model, T; p0 = nothing, v0 = nothing)

Calculate the Characteristic Isobaric Inflection Curve (CIIC) point at a given temperature for a thermodynamic model.

The CIIC represents points of maximum isobaric heat capacity (Cp) along isotherms, in contrast to
Widom lines which represent maxima of Cp along isobars.

# Arguments
- `model`: Thermodynamic model used for calculations.
- `T`: Temperature at which to calculate the CIIC pressure `[K]`.
- `p0`: Optional initial pressure guess `[Pa]`. If not provided, a default value will be used.
- `v0`: Optional initial volume guess `[m³]`. If not provided, a default value will be used.

# Returns
- `p_ciic`: Pressure at the CIIC point for the specified temperature `[Pa]`.
- `v_ciic`: Volume at the CIIC point `[m³]`.

# Example
```julia
model = PCSAFT(["methane"])
T = 200.0  # 200 K
p_ciic, v_ciic = ciic_pressure(model, T)
```
"""
function ciic_pressure(model,T;p0 = nothing,v0 = nothing)
    if !isnothing(p0) && isnothing(v0)
        _v0 = volume(model,p0,T,phase = :l)
    elseif isnothing(p0) && !isnothing(v0)
        _v0 = v0*one(Base.promote_eltype(model,T))
    elseif isnothing(p0) && isnothing(v0)
        crit = crit_pure(model)
        p0 = critical_psat_extrapolation(model,T,crit)
        return ciic_pressure(model,T;p0)
    else
        throw(ArgumentError("cannot specify `v0` and `p0` simultaneously"))
    end
    #=
    ∂Cₚ∂V*∂p∂V: J/K/m3 * Pa/m3
    RTp/(T*v*v)| p = RT/v
    (RT/v)^2/(T*v)
    =#
    lb_v = lb_volume(model,T,SA[1.0])
    scale = abs2(Rgas(model)*T/lb_v)/(T*lb_v)
    function f_ciic(logv)
        v = exp(logv)
        Cₚ(∂v) = VT_isobaric_heat_capacity(model,∂v,T,SA[1.0])
        p(∂v) = pressure(model,∂v,T,SA[1.0])
        ∂Cₚ∂V = Solvers.derivative(Cₚ,v)
        ∂p∂V = Solvers.derivative(p,v)
        return ∂Cₚ∂V*∂p∂V/scale
    end

    v_ciic = exp(Solvers.nlsolve(f_ciic,log(_v0),Roots.Order0()))
    p_ciic = pressure(model,v_ciic,T,SA[1.0])
    return p_ciic,v_ciic
end

"""
    T_ciic, v_ciic = ciic_temperature(model, p; T0 = nothing, v0 = nothing)

Calculate the Characteristic Isobaric Inflection Curve (CIIC) point at a given pressure for a thermodynamic model.

The CIIC represents points of maximum isobaric heat capacity (Cp) along isotherms, in contrast to
Widom lines which represent maxima of Cp along isobars.

# Arguments
- `model`: Thermodynamic model used for calculations.
- `p`: Pressure at which to calculate the CIIC temperature `[Pa]`.
- `T0`: Optional initial temperature guess `[K]`. If not provided, a default value will be used.
- `v0`: Optional initial volume guess `[m³]`. If not provided, a default value will be used.

# Returns
- `T_ciic`: Temperature at the CIIC point for the specified pressure `[K]`.
- `v_ciic`: Volume at the CIIC point `[m³]`.

# Example
```julia
model = PCSAFT(["methane"])
p = 150e5  # 150 bar
T_ciic, v_ciic = ciic_temperature(model, p)
```
"""
function ciic_temperature(model,p;T0 = nothing,v0 = nothing)
    if !isnothing(T0) && isnothing(v0)
        _v0 = volume(model,p,T0,phase = :l)
        _T0 = T0*one(Base.promote_eltype(model,p))
    elseif isnothing(T0) && !isnothing(v0)
        _v0 = v0*one(Base.promote_eltype(model,p))
        crit = crit_pure(model)
        _T0 = critical_tsat_extrapolation(model,p,crit)
    elseif isnothing(T0) && isnothing(v0)
        crit = crit_pure(model)
        T0 = critical_tsat_extrapolation(model,p,crit)
        return ciic_temperature(model,p;T0 = T0)
    else
        _1 = Base.promote_eltype(model,p,T0,v0)
        _v0 = v0*_1
        _T0 = T0*_1
    end

    lb_v = lb_volume(model,_T0,SA[1.0])
    scale = p*p/_T0/lb_volume(model,_T0,SA[1.0])
    #P*P/(T*v)
    function f_ciic(x::AbstractVector)
        v,T = exp(x[1]),x[2]
        Cₚ(∂v) = VT_isobaric_heat_capacity(model,∂v,T,SA[1.0])
        fp(∂v) = pressure(model,∂v,T,SA[1.0])
        ∂Cₚ∂V = Solvers.derivative(Cₚ,v)
        px,∂p∂V = Solvers.f∂f(fp,v)
        ciic_condition = ∂Cₚ∂V*∂p∂V/scale #J*Pa/K/m3/m3
        pressure_eq = (px - p)/p
        return SVector(ciic_condition,pressure_eq)
    end

    #we do a refinement of the initial value first. it improves convergence near the critical point
    function f_ciic(T)
        v = volume(model,p,T,phase = :l,vol0 = _v0)
        return f_ciic(SVector(log(v),T))[1]
    end
    prob = Roots.ZeroProblem(f_ciic,_T0)
    _T00 = Roots.solve(prob,Roots.Order2(),atol = 1e-1)
    _v00 = volume(model,p,_T00,phase = :l, vol0 = _v0)

    res = Solvers.nlsolve2(f_ciic,SVector(log(_v00),_T00),Solvers.Newton2Var())
    v_ciic,T_ciic = exp(res[1]),res[2]
    return T_ciic,v_ciic
end

export widom_pressure, widom_temperature
export ciic_pressure, ciic_temperature
