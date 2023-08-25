module ClapeyronCoolPropExt

using Clapeyron
using CoolProp

function Clapeyron.coolprop_handler()
    CoolProp.CoolProp_jll.libcoolprop
end

#TODO: CoolProp-Clapeyron integration?

end