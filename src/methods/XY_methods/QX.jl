function QT_property(model,q,T,z,f::F,p0) where F
    if f == temperature
        return T
    end

    res = qt_flash(model,q,T,z,p0 = p0)
    if isone(numphases(res)) && !isone(q) !izero(q)
        #What to do here?
    end
    if f == temperature
        return temperature(res)
    elseif f == pressure
        return pressure(res)
    else
        return f(model,res)
    end
end

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
    if f == pressure
        return p
    end

    res = qp_flash(model,q,p,z,T0 = T0)
    if isone(numphases(res)) && !isone(q) !iszero(q)
        #What to do here?
    end
    if f == temperature
        return temperature(res)
    elseif p == pressure
        return pressure(res)
    else
        return f(model,res)
    end
end

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