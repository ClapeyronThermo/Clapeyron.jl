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

The Widom point represents the temperature where maxima of the isobaric heat capacity occurs along an isobar in the supercritical region.
In particular, `widom_pressure` calculates the widom point at a specified temperature.

# Arguments
- `model`: Thermodynamic model used for calculations
- `T`: Temperature at which to calculate the Widom pressure
- `p0`: Optional initial pressure guess. If not provided, a default value will be used
- `v0`: Optional initial volume guess. If not provided, a default value will be used

# Returns
- `pw`: Widom pressure at the specified temperature
- `vw`: Volume at the Widom pressure and specified temperature

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
    
    function f_wp(v)
        Cₚ(∂v,∂T) = VT_isobaric_heat_capacity(model,∂v,∂T,SA[1.0])
        p(∂v,∂T) = pressure(model,∂v,∂T,SA[1.0])
        ∂Cₚ∂V,∂Cₚ∂T = Solvers.gradient2(Cₚ,v,T)
        ∂p∂V,∂p∂T = Solvers.gradient2(p,v,T)
        return ∂Cₚ∂T - ∂Cₚ∂V * ∂p∂T / ∂p∂V
    end

    v_widom = Solvers.nlsolve(f_wp,_v0,Roots.Order0())
    p_widom = pressure(model,v_widom,T,SA[1.0])
    return p_widom,v_widom
end

"""
    Tw, vw = widom_temperature(model, p; T0 = nothing, v0 = nothing)

Calculate the Widom temperature and corresponding volume at a given pressure for a thermodynamic model.

The Widom point represents the temperature where maxima of the isobaric heat capacity occurs along an isobar in the supercritical region.
In particular, `widom_temperature` calculates the widom point at a specified temperature.

# Arguments
- `model`: Thermodynamic model used for calculations
- `p`: Pressure at which to calculate the Widom temperature
- `T0`: Optional initial temperature guess. If not provided, a default value will be used
- `v0`: Optional initial volume guess. If not provided, a default value will be used

# Returns
- `Tw`: Widom temperature at the specified pressure
- `vw`: Volume at the Widom temperature and specified pressure

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
        throw(ArgumentError("cannot specify `v0` and `T0` simultaneously"))
    end
    
    function f_wp(x)
        v,T = x[1],x[2]
        Cₚ(∂v,∂T) = VT_isobaric_heat_capacity(model,∂v,∂T,SA[1.0])
        fp(∂v,∂T) = pressure(model,∂v,∂T,SA[1.0])
        ∂Cₚ∂V,∂Cₚ∂T = Solvers.gradient2(Cₚ,v,T)
        px,dpx = Solvers.fgradf2(fp,v,T)
        ∂p∂V,∂p∂T = dpx
        widom_condition = ∂Cₚ∂T - ∂Cₚ∂V * ∂p∂T / ∂p∂V
        pressure_eq = (px - p)/p
        return SVector(widom_condition,pressure_eq)
    end

    res = Solvers.nlsolve2(f_wp,SVector(_v0,_T0),Solvers.Newton2Var())
    v_widom,T_widom = res
    return T_widom,v_widom
end

"""
    p_ciic, v_ciic = ciic_pressure(model, T; p0 = nothing, v0 = nothing)

Calculate the Characteristic Isobaric Inflection Curve (CIIC) point at a given temperature for a thermodynamic model.

The CIIC represents points of maximum isobaric heat capacity (Cp) along isotherms, in contrast to 
Widom lines which represent maxima of Cp along isobars. 

# Arguments
- `model`: Thermodynamic model used for calculations
- `T`: Temperature at which to calculate the CIIC pressure
- `p0`: Optional initial pressure guess. If not provided, a default value will be used
- `v0`: Optional initial volume guess. If not provided, a default value will be used

# Returns
- `p_ciic`: Pressure at the CIIC point for the specified temperature
- `v_ciic`: Volume at the CIIC point

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
        return widom_pressure(model,T;p0 = p0)
    else
        throw(ArgumentError("cannot specify `v0` and `p0` simultaneously"))
    end
    
    function f_ciic(v)
        Cₚ(∂v,T) = VT_isobaric_heat_capacity(model,∂v,T,SA[1.0])
        p(∂v,T) = pressure(model,∂v,T,SA[1.0])
        ∂Cₚ∂V = Solvers.derivative(Cₚ,v)
        ∂p∂V = Solvers.derivative(p,v)
        return ∂Cₚ∂V*∂p∂V
    end

    v_ciic = Solvers.nlsolve(f_ciic,_v0,Roots.Order0())
    p_ciic = pressure(model,v_widom,T,SA[1.0])
    return p_ciic,v_ciic
end

"""
    T_ciic, v_ciic = ciic_temperature(model, p; T0 = nothing, v0 = nothing)

Calculate the Characteristic Isobaric Inflection Curve (CIIC) point at a given pressure for a thermodynamic model.

The CIIC represents points of maximum isobaric heat capacity (Cp) along isotherms, in contrast to 
Widom lines which represent maxima of Cp along isobars. 

# Arguments
- `model`: Thermodynamic model used for calculations
- `p`: Pressure at which to calculate the CIIC temperature
- `T0`: Optional initial temperature guess. If not provided, a default value will be used
- `v0`: Optional initial volume guess. If not provided, a default value will be used

# Returns
- `T_ciic`: Temperature at the CIIC point for the specified pressure
- `v_ciic`: Volume at the CIIC point

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
        return widom_temperature(model,p;T0 = T0)
    else
        throw(ArgumentError("cannot specify `v0` and `T0` simultaneously"))
    end
    
    function f_ciic(x)
        v,T = x[1],x[2]
        Cₚ(∂v) = VT_isobaric_heat_capacity(model,∂v,T,SA[1.0])
        fp(∂v) = pressure(model,∂v,T,SA[1.0])
        ∂Cₚ∂V = Solvers.derivative(Cₚ,v,T)
        px,∂p∂V = Solvers.d∂f(fp,v)
        ciic_condition = ∂Cₚ∂V*∂p∂V
        pressure_eq = (px - p)/p
        return SVector(ciic_condition,pressure_eq)
    end

    res = Solvers.nlsolve2(f_ciic,SVector(_v0,_T0),Solvers.Newton2Var())
    v_ciic,T_ciic = res
    return T_ciic,v_ciic
end

export widom_pressure, widom_temperature
export ciic_pressure, ciic_temperature
