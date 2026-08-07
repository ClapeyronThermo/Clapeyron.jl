
function VT_property(model,V,T,z,f::F,phase,p0) where {F}
    XX = Base.promote_eltype(model,V,T,z)
    z isa Number && return VT_property(model,V,T,SVector(z),f,phase,p0)
    f == volume && return XX(V)
    f == temperature && return XX(T)
    !is_unknown(phase) && return spec_to_vt(model,V,T,z,f)
    res = vt_flash(model,V,T,z,p0 = p0)
    return f(model,res)
end

module VT
import Clapeyron
for f in Clapeyron.CLAPEYRON_PROPS
    @eval begin
        function $f(model,V,T,z = Clapeyron.SA[1.0],p0 = nothing,phase = :unknown)
            Clapeyron.VT_property(model,V,T,z,Clapeyron.$f,phase,p0)
        end
    end
end
function flash(model,V,T,z = Clapeyron.SA[1.0],args...;kwargs...)
    return Clapeyron.vt_flash(model,V,T,z,args...;kwargs...)
end
end #VT module
