module ClapeyronCoolPropExt

using Clapeyron
using CoolProp

function Clapeyron.coolprop_handler()
    Base.Libc.Libdl.dlopen(CoolProp.CoolProp_jll.libcoolprop;throw_error = false)
end

#TODO: CoolProp-Clapeyron integration?

end