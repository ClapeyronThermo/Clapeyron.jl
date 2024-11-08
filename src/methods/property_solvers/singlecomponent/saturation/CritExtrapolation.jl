"""
    CritExtrapolationSaturation <: SaturationMethod
    CritExtrapolationSaturation(;crit = nothing)

Saturation method that extrapolates the liquid and vapor volumes from the critical point. the extrapolation is defined as:

```
ρₓ = ρc ± sqrt(6* Tc * (Tc - T)/Tc * ∂²p∂ρ∂T / ∂³p∂ρ³)
```

This method is specially useful at calculating the liquid and vapor volumes when T > 0.999Tc.

## References

1. Bell, I. H., & Deiters, U. K. (2023). Superancillary equations for nonpolar pure fluids modeled with the PC-SAFT equation of state. Industrial & Engineering Chemistry Research. [doi:10.1021/acs.iecr.2c02916](https://doi.org/10.1021/acs.iecr.2c02916)
"""
struct CritExtrapolationSaturation{T} <: SaturationMethod
    crit::NTuple{3,T}
end

function CritExtrapolationSaturation(;crit = nothing)
    if crit == nothing
        return CritExtrapolationSaturation{Float64}((0.0,0.0,0.0))
    else
        _crit = promote(crit...)
        T = typeof(_crit[1])
        return CritExtrapolationSaturation{T}(_crit)
    end
end

function saturation_pressure_impl(model::EoSModel, T, method::CritExtrapolationSaturation{C}) where C
    crit = method.crit
    T_c, p_c, V_c = crit
    if T_c == p_c == V_c == 0.0
        return saturation_pressure_impl(model,T,CritExtrapolationSaturation(crit_pure(model)))
    end
    vl,vv = critical_vsat_extrapolation(model,T,T_c,V_c)
    p = pressure(model,vl,T)
    return (p,vl,vv)
end

function saturation_temperature_impl(model::EoSModel, p, method::CritExtrapolationSaturation{C}) where C
    crit = method.crit
    T_c, p_c, V_c = crit
    if T_c == p_c == V_c == 0.0
        return saturation_temperature_impl(model,T,CritExtrapolationSaturation(crit_pure(model)))
    end
    T = critical_tsat_extrapolation(model,p,crit)
    vl,vv = critical_vsat_extrapolation(model,T,crit)
    return T,vl,vv
end

export CritExtrapolationSaturation
