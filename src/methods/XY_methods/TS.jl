function TS_property(model,T,s,z,f::F,phase,p0,threaded) where F
    z isa Number && return TS_property(model,T,s,SVector(z),f,phase,T0,threaded)
    XX = Base.promote_eltype(model,T,s,z)
    f == entropy && return XX(s)
    f == temperature && return XX(T)

    if f == pressure && length(model) == 1
        z1 = SVector(z[1])
        return Pproperty(model,T,s,z1,entropy,p0 = p0,phase = phase,threaded = threaded)
    end

    if !is_unknown(phase)
        p,calc_phase = _Pproperty(model,T,s,z,p0 = p0,phase = phase,threaded = threaded)
        if calc_phase != :eq && calc_phase != :failure
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

module TS
import Clapeyron
for f in Clapeyron.CLAPEYRON_PROPS
    @eval begin
        function $f(model,T,s,z = Clapeyron.SA[1.0];phase = :unknown,p0 = nothing, threaded = true)
            Clapeyron.TS_property(model,T,s,z,Clapeyron.$f,phase,p0,threaded)
        end
    end
end
function flash(model,T,s,z = Clapeyron.SA[1.0],args...;kwargs...)
    return Clapeyron.ts_flash(model,T,s,z,args...;kwargs...)
end
end #module
