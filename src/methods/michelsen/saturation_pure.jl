import Optim: optimize
import Roots: find_zero
using NLsolve: nlsolve, only_fj!

# objective functions for Pmin and Pmax initiation method
function fobj_pmax(model::EoSModel, V, T, z=[1.])
    return - pressure(model, V, T, z)
end

function fobj_pmin(model::EoSModel, V, T, z=[1.])
    return pressure(model, V, T, z)
end

function psat_init(model::EoSModel, T, Tc, Vc)
    # Function to get an initial guess for the saturation pressure at a given temperature
    z = [1.]
    RT = R̄*T
    Tr = T/Tc
    # Zero pressure initiation
    if Tr < 0.8
        P0 = 0.
        vol_liq0 = volume(model, P0, T, phase=:liquid)
        ares = a_res(model, vol_liq0, T, z)
        lnϕ_liq0 = ares - 1. + log(RT/vol_liq0)
        P0 = exp(lnϕ_liq0)
    # Pmin, Pmax initiation
    elseif Tr <= 1.0
        low_v = Vc
        up_v = 5 * Vc
        fmax(V) = fobj_pmax(model, V, T)
        sol_max = optimize(fmax, low_v, up_v)
        v_max = sol_max.minimizer
        P_max = - sol_max.minimum

        low_v = 1e-10*Vc
        up_v = Vc
        fmin(V) = fobj_pmin(model, V, T)
        sol_min = optimize(fmin, low_v, up_v)
        v_min = sol_min.minimizer
        P_min = sol_min.minimum
        P0 = (max(0., P_min) + P_max) / 2
    else
        P0 = NaN
    end
    return  P0
end


function psat_fugacity(model::EoSModel, T, p0, vol0=[nothing, nothing])
    # Objetive function to solve saturation pressure using the pressure as iterable variable
    # T = Saturation Temperature
    # p0 = initial guess for the saturation pressure
    # vol0 = initial guesses for the phase volumes = [vol liquid, vol vapor]
    # out = Saturation Pressure, vol liquid, vol vapor
    vol_liq0, vol_vap0 = vol0

    RT = R̄*T
    P = 1. * p0
    # Solving the phase volumes for the first iteration
    vol_liq = volume(model, P, T, phase=:liquid, vol0=vol_liq0)
    vol_vap = volume(model, P, T, phase=:vapor, vol0=vol_vap0)

    itmax = 20
    for i in 1:itmax
        # Computing chemical potential
        μ_liq = VT_chemical_potential_res(model, vol_liq, T)[1]
        μ_vap = VT_chemical_potential_res(model, vol_vap, T)[1]

        Z_liq = P*vol_liq/RT
        Z_vap = P*vol_vap/RT

        lnϕ_liq = μ_liq/RT - log(Z_liq)
        lnϕ_vap = μ_vap/RT - log(Z_vap)
        # Updating the saturation pressure
        FO = lnϕ_vap - lnϕ_liq
        dFO = (Z_vap - Z_liq) / P
        dP = FO / dFO
        P = P - dP
        if abs(dP) < 1e-8; break; end
        # Updating the phase volumes
        vol_liq = volume(model, P, T, phase=:liquid, vol0=vol_liq)
        vol_vap = volume(model, P, T, phase=:vapor, vol0=vol_vap)
    end
    return P, vol_liq, vol_vap
end


function ∂Helmholtz(model::EoSModel, ρ, T, z=[1.0])
    # Auxiliar functions to compute Helmholtz Energy and its first density derivative
    f(dρ) = eos(model, 1. / dρ, T, z)
    A, ∂A, = Solvers.f∂f(f, ρ)
    return A, ∂A
end

function ∂2Helmholtz(model::EoSModel, ρ, T, z=[1.0])
    # Auxiliar functions to compute Helmholtz Energy and its first and second density derivatives
    f(dρ) = eos(model, 1. / dρ, T, z)
    A, ∂A, ∂2A = Solvers.f∂f∂2f(f, ρ)
    return A, ∂A, ∂2A
end


function fobj_psat!(model::EoSModel, ρ, T, F, J)
    # Objetive function to solve saturation pressure using the densities as iterable variables
    # T = Saturation Temperature
    # ρ = initial guess for the phases densities = [ρ liquid, ρ vapor]
    # F = vector for objective function
    # J = matriz for objective function jacobian

    ρ_liq, ρ_vap = ρ

    if !(J == nothing)
        A_liq, ∂A_liq, ∂2A_liq = ∂2Helmholtz(model, ρ_liq, T)
        ∂P_liq = 2. * ρ_liq * ∂A_liq + ρ_liq^2 * ∂2A_liq
        ∂μ_liq = ρ_liq*∂2A_liq + 2*∂A_liq

        A_vap, ∂A_vap, ∂2A_vap = ∂2Helmholtz(model, ρ_vap, T)
        ∂P_vap = 2. * ρ_vap * ∂A_vap + ρ_vap^2 * ∂2A_vap
        ∂μ_vap = ρ_vap*∂2A_vap + 2*∂A_vap
        J[1, 1] = ∂μ_liq
        J[1, 2] = -∂μ_vap
        J[2, 1] = ∂P_liq
        J[2, 2] = -∂P_vap
    else
        A_liq, ∂A_liq = ∂Helmholtz(model, ρ_liq, T)
        A_vap, ∂A_vap = ∂Helmholtz(model, ρ_vap, T)
    end

    P_liq = ∂A_liq*ρ_liq^2
    μ_liq = A_liq + ρ_liq * ∂A_liq

    P_vap = ∂A_vap*ρ_vap^2
    μ_vap = A_vap + ρ_vap * ∂A_vap

    if !(F == nothing)
        F[1] = μ_liq - μ_vap
        F[2] = P_liq - P_vap
    end

end


function psat(model::EoSModel, T; p0=nothing, vol0=[nothing, nothing])
    # Function to solve saturation pressure of a pure fluid
    # T = Saturation Temperature
    # p0 = initial guess for the saturation pressure
    # vol0 = initial guesses for the phase volumes = [vol liquid, vol vapor]
    # out = Saturation Pressure, vol liquid, vol vapor

    init = false

    vol_liq0, vol_vap0 = vol0

    if p0 == nothing && vol_liq0 != nothing && vol_vap0 != nothing
        # if no initial pressure is given solve by using equality of chemical potentials and pressure
        method = :chempot
    elseif p0 == nothing && vol_liq0 == nothing && vol_vap0 == nothing
        # if no initial guess is given estimate by using psat_init function
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
        p0 = psat_init(model, T, Tc, Vc)
    end

    if method == :fug
        P, vol_liq, vol_vap = psat_fugacity(model, T, p0, vol0)
    elseif method == :chempot
        ρ0 = [1/vol_liq0, 1/vol_vap0]
        ofpsat(F, J, ρ) = fobj_psat!(model, ρ, T, F, J)
        sol = nlsolve(only_fj!(ofpsat), ρ0, method = :newton)
        ρ = sol.zero
        vol_liq, vol_vap = 1 ./ ρ
        P = pressure(model, vol_vap, T, [1.])
    end

    return P, vol_liq, vol_vap
end

############ Solving saturation pressure at given pressure
function obj_tsat(model::EoSModel, T, P)
    global vol_liq
    global vol_vap

    RT = R̄*T

    vol_liq = volume(model, P, T, phase=:liquid, vol0=vol_liq)
    vol_vap = volume(model, P, T, phase=:vapor, vol0=vol_vap)

    μ_liq = VT_chemical_potential_res(model, vol_liq, T)[1]
    μ_vap = VT_chemical_potential_res(model, vol_vap, T)[1]

    Z_liq = P*vol_liq/RT
    Z_vap = P*vol_vap/RT

    lnϕ_liq = μ_liq/RT - log(Z_liq)
    lnϕ_vap = μ_vap/RT - log(Z_vap)
    FO = lnϕ_vap - lnϕ_liq
    return FO
end

function tsat(model::EoSModel, P, T0)
    # Function to solve saturation temperature of a pure fluid
    # P = Saturation Pressure
    # T0 = initial guess for the saturation temperature
    # out = Saturation Temperature, vol liquid, vol vapor
    global vol_liq
    global vol_vap

    vol_liq = nothing
    vol_vap = nothing

    ftsat(T) = obj_tsat(model, T, P)
    T = find_zero(ftsat, T0)
    return T, vol_liq, vol_vap
end
