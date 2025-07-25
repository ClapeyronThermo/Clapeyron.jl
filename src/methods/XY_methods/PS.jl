function PS_property(model,p,s,z,f::F,phase,T0,threaded) where F
    z isa Number && return PS_property(model,p,s,SVector(z),f,phase,T0,threaded)
    XX = Base.promote_eltype(model,p,s,z)
    f == entropy && return XX(s)
    f == pressure && return XX(p)

    if f == temperature && length(model) == 1
        z1 = SVector(z[1])
        return Tproperty(model,p,s,z1,entropy,T0 = T0,phase = phase,threaded = threaded)
    end

    if !is_unknown(phase)
        T,calc_phase = _Tproperty(model,p,h,z,entropy,T0 = T0,phase = phase,threaded = threaded)
        if calc_phase != :eq && calc_phase != :failure
            return f(model,p,T,z;phase = calc_phase)
        elseif calc_phase == :eq
            supports_lever_rule(f) || thow(invalid_property_multiphase_error(f))
            result = ps_flash(model,p,s,z,T0)
            return f(model,result)
        else
            return f(model,p,T,z;phase = phase)
        end
    end
    res = ps_flash(model,p,s,z,T0 = T0)
    return f(model,res)
end

module PS
import Clapeyron
for f in Clapeyron.CLAPEYRON_PROPS
    @eval begin
        function $f(model,p,s,z = Clapeyron.SA[1.0];phase = :unknown,T0 = nothing, threaded = true)
            Clapeyron.PS_property(model,p,s,z,Clapeyron.$f,phase,T0,threaded)
        end
    end
end
function flash(model,p,s,z = Clapeyron.SA[1.0],args...;kwargs...)
    return Clapeyron.ps_flash(model,p,s,z,args...;kwargs...)
end
end #module