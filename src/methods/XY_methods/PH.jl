function PH_property(model,p,h,z,f::F,phase,T0,threaded) where F
    XX = Base.promote_eltype(model,p,h,z)
    if f == enthalpy
        return XX(h)
    end

    if f == pressure
        return XX(p)
    end

    if !is_unknown(phase)
        T,calc_phase = _Tproperty(model,p,h,z,T0 = T0,phase = phase,threaded = threaded)
        if calc_phase != :eq && calc_phase != :failure
            return XX(f(model,p,T,z;phase = calc_phase))
        elseif calc_phase == :eq && !supports_lever_rule(f)
            thow(invalid_property_multiphase_error(f))
        elseif calc_phase == :eq && supports_lever_rule(f)
            result = ph_flash(model,p,h,z,T0)
            return XX(f(model,result))
        else
            return XX(f(model,p,T,z;phase = phase))
        end
    else
        res = ph_flash(model,p,h,z,T0 = T0)
        if f == temperature
            return XX(temperature(res))
        else
            return XX(f(model,res))
        end
    end
end

module PH
import Clapeyron
for f in Clapeyron.CLAPEYRON_PROPS
    @eval begin
        function $f(model,p,h,z = Clapeyron.SA[1.0];phase = :unknown,T0 = nothing, threaded = true)
            Clapeyron.PH_property(model,p,h,z,Clapeyron.$f,phase,T0,threaded)
        end
    end
end
function flash(model,p,h,z = Clapeyron.SA[1.0],args...;kwargs...)
    return Clapeyron.ph_flash(model,p,h,z,args...;kwargs...)
end
end  #module