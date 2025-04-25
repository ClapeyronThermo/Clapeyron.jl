
function VT_property(model,V,T,z,f::F,phase,p0) where {F}
    if z isa Number
        return VT_property(model,V,T,SA[z],f,phase,p0)
    end

    if f == volume
        return V
    end

    if f == temperature
        return T
    end

    if !is_unknown(phase)
        return spec_to_vt(model,V,T,z,f)
    end

    res = vt_flash(model,V,T,z,p0 = p0)
    if f == temperature
        return temperature(res)
    elseif p == pressure
        return pressure(res)
    else
        return f(model,res)
    end
end



module VT
import Clapeyron
for f in CLAPEYRON_PROPS
    @eval begin
        function $f(model,V,T,z = Clapeyron.SA[1.0],p0 = nothing,phase = :unknown)
            Clapeyron.VT_property(model,V,T,z,Clapeyron.$f,p0,phase)
        end

        function flash(model,V,T,z = Clapeyron.SA[1.0],args...;kwargs...)
            return Clapeyron.vt_flash(model,V,T,z,args...;kwargs...)
        end
    end
end
end #VT module
