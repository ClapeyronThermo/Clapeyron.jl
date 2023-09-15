module ClapeyronCoolPropExt

using Clapeyron
using CoolProp

function Clapeyron.coolprop_handler()
    Base.Libc.Libdl.dlopen(CoolProp.CoolProp_jll.libcoolprop;throw_error = false)
end

#TODO: CoolProp-Clapeyron integration?

"""
    coolprop_crit_data(components)
returns a named tuple with critical and molecular weight data extracted from CoolProp.
"""
function Clapeyron.coolprop_crit_data(components)
    comps = Clapeyron.format_components(components)
    Tc = zeros(length(comps))
    Vc = zeros(length(comps))
    acentricfactor = zeros(length(comps))
    Pc = zeros(length(comps))
    Mw = zeros(length(comps))

    for (i,comp) in pairs(comps)
        Tc[i] = PropsSI("Tcrit",comp)
        Pc[i] = PropsSI("pcrit",comp)
        Vc[i] = 1/PropsSI("rhocrit",comp)
        Mw[i] = 1000*PropsSI("molarmass",comp)
        acentricfactor[i] =PropsSI("acentric",comp)
    end
    return (;Mw,Tc,Pc,Vc,acentricfactor)

end


end