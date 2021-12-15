"""
    acentric_factor(model::EoSModel)

calculates the acentric factor using its definition:

    ω = -log10(psatᵣ) -1, at Tᵣ = 0.7
To do so, it calculates the critical temperature (using `crit_pure`) and performs a saturation calculation (with `sat_pure`)

"""
function acentric_factor(model)
    T_c,p_c,_ = crit_pure(model)
    T = 0.7*T_c
    p = first(saturation_pressure(model,T))
    p_r = p/p_c
    return -log10(p_r) - 1.0
end

include("saturation.jl")
include("crit_pure.jl")

export saturation_pressure, saturation_temperature, crit_pure, enthalpy_vap