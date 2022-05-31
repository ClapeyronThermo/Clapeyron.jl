
function psat_init(model::EoSModel, T, Tc, Vc)
    # Function to get an initial guess for the saturation pressure at a given temperature
    z = SA[1.] #static vector
    _0 = zero(T+Tc+Vc)
    RT = R̄*T
    Tr = T/Tc
    # Zero pressure initiation
    if Tr < 0.8
        P0 = _0
        vol_liq0 = volume(model, P0, T, phase=:liquid)
        ares = a_res(model, vol_liq0, T, z)
        lnϕ_liq0 = ares - 1. + log(RT/vol_liq0)
        P0 = exp(lnϕ_liq0)
    # Pmin, Pmax initiation
    elseif Tr <= 1.0
        low_v = Vc
        up_v = 5 * Vc
        #note: P_max is the pressure at the maximum volume, not the maximum pressure
        fmax(V) = pressure(model, V, T)
        sol_max = Solvers.optimize(fmax, (low_v, up_v))
        v_max = Solvers.x_sol(sol_max)
        P_max = Solvers.x_minimum(sol_max)

        low_v = lb_volume(model)
        up_v = Vc

        #note: P_min is the pressure at the minimum volume, not the minimum pressure
        fmin(V) = -pressure(model, V, T)
        sol_min = Solvers.optimize(fmin, (low_v,up_v))
        v_min = Solvers.x_sol(sol_min)
        P_min = Solvers.x_minimum(sol_min)
        P0 = (max(0., P_min) + P_max) / 2
    else
        P0 = _0/_0 #NaN, but propagates the type
    end
    return  P0
end

struct IsoFugacitySaturation{T} <: SaturationMethod
    p0::T
    vl::Union{Nothing,T}
    vv::Union{Nothing,T}
    max_iters::Int
    p_tol::Float64
end

function IsoFugacitySaturation(;p0 = nothing,vl = nothing,vv = nothing,max_iters = 20,p_tol = sqrt(eps(Float64)))
    p0 === nothing && (p0 = NaN)
    if vl !== nothing
        p0,vl = promote(p0,vl)
    elseif vv !== nothing
        p0,vv = promote(p0,vv)
    elseif (vv !== nothing) & (vl !== nothing)
        p0,vl,vv = promote(p0,vl,vv)
    else
    end
    return IsoFugacitySaturation(p0,vl,vv,max_iters,p_tol)
end

function saturation_pressure_impl(model::EoSModel,T,method::IsoFugacitySaturation)
    vol0 = (method.vl,method.vv,T)
    p0 = method.p0
    if isnan(p0)
        Tc, Pc, Vc = crit_pure(model)
        if Tc < T
            nan = p0*zero(T)/zero(T)
            return (nan,nan,nan)
        end
        p0 = psat_init(model, T, Tc, Vc)
    end
    return psat_fugacity(model,T,p0,vol0,method.max_iters,method.p_tol)
end

function psat_fugacity(model::EoSModel, T, p0, vol0=(nothing, nothing),max_iters = 20,p_tol = sqrt(eps(Float64)))
    # Objetive function to solve saturation pressure using the pressure as iterable variable
    # T = Saturation Temperature
    # p0 = initial guess for the saturation pressure
    # vol0 = initial guesses for the phase volumes = [vol liquid, vol vapor]
    # out = Saturation Pressure, vol liquid, vol vapor
    z = SA[1.]
    RT = R̄*T
    P = 1. * p0
    vol_liq0, vol_vap0 = vol0
    vol_liq0 === nothing && (vol_liq0 = x0_volume_liquid(model,T,z))
    vol_vap0 === nothing && (vol_vap0 = x0_volume_gas(model,P,T,z))

    # Solving the phase volumes for the first iteration

    vol_liq = _volume_compress(model, P, T, z, vol_liq0)
    vol_vap = _volume_compress(model, P, T, z, vol_vap0)

    itmax = max_iters
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
        if abs(dP) < p_tol; break; end
        # Updating the phase volumes
        vol_liq = _volume_compress(model, P, T, z,vol_liq)
        vol_vap = _volume_compress(model, P, T, z,vol_vap)
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

    if !(J === nothing)
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

    if !(F === nothing)
        F[1] = μ_liq - μ_vap
        F[2] = P_liq - P_vap
    end
end

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
        return psat_fugacity(model, T, p0, vol0)
    elseif method == :chempot
        return psat_chempot(model,T,vol_liq0,vol_vap0)
    end
    return P, vol_liq, vol_vap
end

struct ChemPotDensitySaturation{T} <: SaturationMethod
    vl::Union{Nothing,T}
    vv::Union{Nothing,T}
end

function ChemPotDensitySaturation(;vl = nothing,vv = nothing)
    if (vl === nothing) && (vv === nothing)
        return ChemPotVSaturation{Nothing}(nothing,nothing)
    elseif !(vl === nothing) && (vv === nothing)
        vl = float(vl)
        return ChemPotVSaturation(vl,vv)
    elseif (vl === nothing) && !(vv === nothing)
        vv = float(vv)
        return ChemPotVSaturation(vl,vv)
    else
        T = one(vl)/one(vv)
        vl,vv,_ = promote(vl,vv,T)
        return ChemPotVSaturation(vl,vv)
    end
end

function saturation_pressure_impl(model::EoSModel, T, method::ChemPotDensitySaturation{Nothing})
    x0 = x0_sat_pure(model,T) .|> exp10
    vl,vv = x0
    return saturation_pressure_impl(model,T,ChemPotVSaturation(vl,vv))
end

function saturation_pressure_impl(model::EoSModel,T,method::ChemPotDensitySaturation)
    return psat_chempot(model,T,method.vl,method.vv)
end

function psat_chempot(model,T,vol_liq0,vol_vap0)
    ρl0 = 1/vol_liq0
    ρv0 = 1/vol_vap0
    ρ0 = vec2(ρl0,ρv0,T)
    ofpsat(F, J, ρ) = fobj_psat!(model, ρ, T, F, J)
    # sol = NLsolve.nlsolve(only_fj!(ofpsat), ρ0, method = :newton)
    sol = Solvers.nlsolve(Solvers.only_fj!(ofpsat), ρ0, LineSearch(Newton()))
    ρ = Solvers.x_sol(sol)
    vol_liq, vol_vap = 1 ./ ρ
    P = pressure(model, vol_vap, T)
    return P, vol_liq, vol_vap
end

############ Solving saturation pressure at given pressure
function obj_tsat(model::EoSModel, T, P)
    global vol_liq
    global vol_vap
    #to avoid using globals, the best thing here is using a cache (TODO)
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

    T = Roots.find_zero(ftsat, T0)
    return T, vol_liq, vol_vap
end
