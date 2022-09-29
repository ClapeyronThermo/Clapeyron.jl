#phasepy api
function psat(model::EoSModel, T; p0=nothing, vol0=(nothing, nothing))
    # Function to solve saturation pressure of a pure fluid
    # T = Saturation Temperature
    # p0 = initial guess for the saturation pressure
    # vol0 = initial guesses for the phase volumes = [vol liquid, vol vapor]
    # out = Saturation Pressure, vol liquid, vol vapor

    T = T*one(T)/one(T)
    init = false

    vol_liq0, vol_vap0 = vol0

    if p0 == nothing && vol_liq0 != nothing && vol_vap0 != nothing
        # if no initial pressure is given solve by using equality of chemical potentials and pressure
        method = :chempot
    elseif p0 == nothing && vol_liq0 == nothing && vol_vap0 == nothing
        # if no initial guess is given estimate by using x0_psat function
        init = true
        method = :fug
    elseif p0 != nothing # && vol_liq0 == nothing && vol_vap0 == nothing
        # if there is an initial guess for fugacity then use isofugacity method
        method = :fug
    end

    if init
        # In the future the critical point should be computed only once
        # and then reused
        Tc, Pc, Vc = crit_pure(model)
        p0 = x0_psat(model, T, Tc, Vc)
    end

    if method == :fug
        return psat_fugacity(model, T, p0, vol0)
    elseif method == :chempot
        return psat_chempot(model,T,vol_liq0,vol_vap0)
    end
    return P, vol_liq, vol_vap
end

############ Solving saturation pressure at given pressure
function tsat(model::EoSModel, P, T0)
    # Function to solve saturation temperature of a pure fluid
    # P = Saturation Pressure
    # T0 = initial guess for the saturation temperature
    # out = Saturation Temperature, vol liquid, vol vapor

    return saturation_temperature(model,P,T0)
end
