function QT_property(model,q,T,z,f::F,p0) where F
    XX = Base.promote_eltype(model,q,T,z)
    f == temperature && return XX(T)
    res = qt_flash(model,q,T,z,p0 = p0)
    if isone(numphases(res)) && !isone(q) !iszero(q)
        #What to do here?
    end
    return f(model,res)
end

"""
    QT

Module that stores Clapeyron properties in (molar) vapour fraction - temperature basis.

All bulk properties have the following form:

```julia
property(model,q,t,z;p0 = nothing)
```

A vapour fraction - temperature flash is done to check if the input pair corresponds to one or more phases.
a `p0` argument can be used to provide an initial pressure guess to the Q-T flash.
"""
module QT
import Clapeyron
for f in Clapeyron.CLAPEYRON_PROPS
    @eval begin
        function $f(model,q,T,z = Clapeyron.SA[1.0],p0 = nothing)
            Clapeyron.QT_property(model,q,T,z,Clapeyron.$f,p0)
        end
    end
end

function flash(model,q,T,z = Clapeyron.SA[1.0],args...;kwargs...)
    return Clapeyron.qt_flash(model,q,T,z,args...;kwargs...)
end
end #module

function QP_property(model,q,p,z,f::F,T0) where F
    XX = Base.promote_eltype(model,q,p,z)
    f == pressure && return XX(p)
    res = qp_flash(model,q,p,z,T0 = T0)
    if isone(numphases(res)) && !isone(q) !iszero(q)
        #What to do here?
    end
    return f(model,res)
end

"""
    QP

Module that stores Clapeyron properties in (molar) vapour fraction - pressure basis.

All bulk properties have the following form:

```julia
property(model,q,p,z;p0 = nothing)
```

A vapour fraction - pressure flash is done to check if the input pair corresponds to one or more phases.
a `T0` argument can be used to provide an initial temperature guess to the Q-P flash.
"""
module QP
import Clapeyron
for f in Clapeyron.CLAPEYRON_PROPS
    @eval begin
        function $f(model,q,p,z = Clapeyron.SA[1.0],T0 = nothing)
            Clapeyron.QP_property(model,q,p,z,Clapeyron.$f,T0)
        end
    end
end
function flash(model,q,p,z = Clapeyron.SA[1.0],args...;kwargs...)
    return Clapeyron.qp_flash(model,q,p,z,args...;kwargs...)
end
end #module