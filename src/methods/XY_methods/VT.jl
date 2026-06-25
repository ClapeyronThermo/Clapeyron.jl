
function VT_property(model,V,T,z,f::F,p0) where {F}
    XX = Base.promote_eltype(model,V,T,z)
    z isa Number && return VT_property(model,V,T,SVector(z),f,p0)
    f == volume && return XX(V)
    f == temperature && return XX(T)
    res = vt_flash(model,V,T,z,p0 = p0)
    return f(model,res)
end

"""
    VT

Module that stores Clapeyron properties in (total) volume-temperature basis.

All bulk properties have the following form:

```julia
property(model,V,T,z;p0 = nothing)
```

A volume-temperature flash is done to check if the volume-temperature input pair corresponds to one or more phases.
a `p0` argument can be used to provide an initial pressure guess to the V-T flash.
To evaluate the property directly in the V-T base, use the `VT0` module instead.
"""
module VT
import Clapeyron
for f in Clapeyron.CLAPEYRON_PROPS
    @eval begin
        function $f(model,V,T,z = Clapeyron.SA[1.0];p0 = nothing)
            Clapeyron.VT_property(model,V,T,z,Clapeyron.$f,p0)
        end
    end
end
function flash(model,V,T,z = Clapeyron.SA[1.0],args...;kwargs...)
    return Clapeyron.vt_flash(model,V,T,z,args...;kwargs...)
end
end #VT module