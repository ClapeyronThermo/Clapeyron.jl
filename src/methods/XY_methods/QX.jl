function QT_property(model, q, T, z, f::F, p0) where F
    XX = Base.promote_eltype(model, q, T, z)
    f == temperature && return XX(T)
    res = qt_flash(model, q, T, z, p0=p0)
    if isone(numphases(res)) && !isone(q)
        !iszero(q)
        #What to do here?
    end
    return f(model, res)
end

"""
    QT

Module that stores Clapeyron properties in (molar) vapour fraction - temperature basis.

All bulk properties have the following form:

```julia
property(model, q, t, z; p0=nothing)
```

A vapour fraction - temperature flash is done to check if the input pair corresponds to one or more phases.
a `p0` argument can be used to provide an initial pressure guess to the Q-T flash.
"""
module QT
import Clapeyron
import Clapeyron: ustrip, uzstrip, with_output_unit, unit_system
for f in Clapeyron.CLAPEYRON_PROPS
    VT_f = Symbol(:VT_, f)
    @eval begin
        function $f(model, q, T, z=Clapeyron.SA[1.0], p0=nothing, output=nothing)
            T̄, z̄, p̄0 = ustrip(T, temperature), uzstrip(model, z), ustrip(p0, pressure)
            prop = Clapeyron.QT_property(model, q, T̄, z̄, Clapeyron.$f, p̄0)
            return with_output_unit(prop, (unit_system(T, p0, z, output), output), Clapeyron.$VT_f)
        end
    end
end

function flash(model, q, T, z=Clapeyron.SA[1.0], args...; kwargs...)
    return Clapeyron.qt_flash(model, q, T, z, args...; kwargs...)
end
end #module

function QP_property(model, q, p, z, f::F, T0) where F
    XX = Base.promote_eltype(model, q, p, z)
    f == pressure && return XX(p)
    res = qp_flash(model, q, p, z, T0=T0)
    if isone(numphases(res)) && !isone(q)
        !iszero(q)
        #What to do here?
    end
    return f(model, res)
end

"""
    QP

Module that stores Clapeyron properties in (molar) vapour fraction - pressure basis.

All bulk properties have the following form:

```julia
property(model, q, p, z; p0=nothing)
```

A vapour fraction - pressure flash is done to check if the input pair corresponds to one or more phases.
a `T0` argument can be used to provide an initial temperature guess to the Q-P flash.
"""
module QP
import Clapeyron
import Clapeyron: ustrip, uzstrip, with_output_unit, unit_system
for f in Clapeyron.CLAPEYRON_PROPS
    VT_f = Symbol(:VT_, f)
    @eval begin
        function $f(model, q, p, z=Clapeyron.SA[1.0], T0=nothing, output=nothing)
            T̄0, z̄, p̄ = ustrip(T0, temperature), uzstrip(model, z), ustrip(p, pressure)
            prop = Clapeyron.QP_property(model, q, p̄, z̄, Clapeyron.$f, T̄0)
            return with_output_unit(prop, (unit_system(T0, p, z, output), output), Clapeyron.$VT_f)
        end
    end
end
function flash(model, q, p, z=Clapeyron.SA[1.0], args...; kwargs...)
    return Clapeyron.qp_flash(model, q, p, z, args...; kwargs...)
end
end #module
