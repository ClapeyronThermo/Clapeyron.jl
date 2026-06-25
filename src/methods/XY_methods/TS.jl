function TS_property(model,T,s,z,f::F,phase,p0) where F
    z isa Number && return TS_property(model,T,s,SVector(z),f,phase,T0)
    XX = Base.promote_eltype(model,T,s,z)
    f == entropy && return XX(s)
    f == temperature && return XX(T)

    if f == pressure && length(z) == 1
        length(model) == 1 || throw(DimensionMismatch("model and composition vector sizes are inconsistent"))
        z1 = SVector(z[1])
        return Pproperty(model,T,s,z1,entropy,p0 = p0,phase = phase,threaded = false)
    end

    if !is_unknown(phase)
        p,calc_phase = _Pproperty(model,T,s,z,p0 = p0,phase = phase,threaded = false)
        if calc_phase != :eq && calc_phase != :failure
            f == pressure && return XX(p)
            return f(model,p,T,z;phase = calc_phase)
        elseif calc_phase == :eq
            supports_lever_rule(f) || thow(invalid_property_multiphase_error(f))
            result = ts_flash(model,T,s,z,p0)
            return f(model,result)
        else
            return f(model,p,T,z;phase = phase)
        end
    end
    res = ts_flash(model,T,s,z,p0 = p0)
    return f(model,res)
end

"""
    TS

Module that stores Clapeyron properties in temperature - (total) entropy basis.

All bulk properties have the following form:

```julia
property(model,t,s,z;phase = :unknown, p0 = nothing)
```

A temperature-entropy flash is done to check if the input pair corresponds to one or more phases.
A `p0` argument can be used to provide an initial pressure guess to the T-S flash.
If a `phase` argument is specified, then it will be used to skip the flash and instead solve for the input conditions instead.
"""
module TS
import Clapeyron
for f in Clapeyron.CLAPEYRON_PROPS
    @eval begin
        function $f(model,T,s,z = Clapeyron.SA[1.0];phase = :unknown,p0 = nothing)
            Clapeyron.TS_property(model,T,s,z,Clapeyron.$f,phase,p0)
        end
    end
end
function flash(model,T,s,z = Clapeyron.SA[1.0],args...;kwargs...)
    return Clapeyron.ts_flash(model,T,s,z,args...;kwargs...)
end
end #module