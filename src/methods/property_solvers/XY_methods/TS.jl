function TS_property(model,T,s,z,f::F,phase,p0,threaded) where F
    if f == entropy
        return s
    end

    if f == temperature
        return T
    end

    if !is_unknown(phase)
        T,calc_phase = _Tproperty(model,T,s,z,p0 = p0,phase = phase,threaded = threaded)
        if calc_phase != :eq && calc_phase != :failure
            return f(model,p,T,z;phase = calc_phase)
        elseif calc_phase == :eq && !supports_lever_rule(f)
            thow(invalid_property_multiphase_error(f))
        elseif calc_phase == :eq && supports_lever_rule(f)
            result = ts_flash(model,T,s,z,p0)
            return f(model,result)
        else
            return f(model,p,T,z;phase = phase)
        end
    else
        res = ts_flash(model,T,s,z,p0 = p0)
        if f == pressure
            return pressure(res)
        else
            return f(model,res)
        end
    end
end

module TS
import Clapeyron
for f in CLAPEYRON_PROPS
    @eval begin
        function $f(model,T,s,z = Clapeyron.SA[1.0];phase = :unknown,p0 = nothing, threaded = true)
            Clapeyron.TS_property(model,T,s,z,Clapeyron.$f,phase,p0,threaded)
        end
    end

    function flash(model,T,s,z = Clapeyron.SA[1.0],args...;kwargs...)
        return Clapeyron.ts_flash(model,T,s,z,args...;kwargs...)
    end
end
end #module
