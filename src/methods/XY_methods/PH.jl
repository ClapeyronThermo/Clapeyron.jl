function PH_property(model,p,h,z,f::F,phase,T0) where F
    z isa Number && return PH_property(model,p,h,SVector(z),f,phase,T0)
    XX = Base.promote_eltype(model,p,h,z)
    f == enthalpy && return XX(h)
    f == pressure && return XX(p)

    if f == temperature && length(z) == 1
        length(model) == 1 || throw(DimensionMismatch("model and composition vector sizes are inconsistent"))
        z1 = SVector(z[1])
        return Tproperty(model,p,h,z1,enthalpy,T0 = T0,phase = phase,threaded = false)
    end

    if !is_unknown(phase)
        T,calc_phase = _Tproperty(model,p,h,z,enthalpy,T0 = T0,phase = phase,threaded = false)
        if calc_phase != :eq && calc_phase != :failure
            f == temperature && return XX(T)
            return f(model,p,T,z;phase = calc_phase)
        elseif calc_phase == :eq
            supports_lever_rule(f) || thow(invalid_property_multiphase_error(f))
            result = ph_flash(model,p,h,z,T)
            return f(model,result)
        else
            return f(model,p,T,z;phase = phase)
        end
    end
    res = ph_flash(model,p,h,z,T0 = T0)
    return f(model,res)
end

"""
    PH

Module that stores Clapeyron properties in pressure - (total) enthalpy basis.

All bulk properties have the following form:

```julia
property(model,p,h,z;phase = :unknown, T0 = nothing)
```

A pressure-enthalpy flash is done to check if the input pair corresponds to one or more phases.
A `T0` argument can be used to provide an initial temperature guess to the P-H flash.
If a `phase` argument is specified, then it will be used to skip the flash and instead solve for the input conditions instead.
"""
module PH
import Clapeyron
for f in Clapeyron.CLAPEYRON_PROPS
    @eval begin
        function $f(model,p,h,z = Clapeyron.SA[1.0];phase = :unknown,T0 = nothing)
            Clapeyron.PH_property(model,p,h,z,Clapeyron.$f,phase,T0)
        end
    end
end
function flash(model,p,h,z = Clapeyron.SA[1.0],args...;kwargs...)
    return Clapeyron.ph_flash(model,p,h,z,args...;kwargs...)
end
end  #module