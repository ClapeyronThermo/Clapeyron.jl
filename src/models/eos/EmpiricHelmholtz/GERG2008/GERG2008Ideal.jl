struct GERG2008Ideal <: IdealModel
    lengthcomponents::Int
    icomponents::Vector{String}
    acentric_factor::Vector{Float64}
    molecularWeight::Vector{Float64}
    criticalDensity::Vector{Float64}
    criticalVolume::Vector{Float64}
    criticalTemperature::Vector{Float64}
    criticalPressure::Vector{Float64}
    ideal_iters::Vector{Vector{Int64}}
    nr::Matrix{Float64}
    zeta::Matrix{Float64}
end



function GERG2008Ideal(components::Vector{String};idealmodel = GER2008Ideal,verbose=false)
lengthcomponents = length(components)
icomponents = copy(components)
xsel = GERG2008_xsel(components)

    return GERG2008Ideal(
        lengthcomponents,
        icomponents,
        acentric_factor,
        molecularWeight,
        criticalDensity,
        criticalVolume,
        criticalTemperature,
        criticalPressure,
        ideal_iters,
        nr,
        zeta
    )
end

function _f0(model::GERG2008Ideal, ρ, T, x)
    if model.lengthcomponents == 1
        return __f0(model,ρ,T)
    else
        return __f0(model,ρ,T,x)
    end
end


function __f0(model::GERG2008Ideal, ρ, T, x)


    RR = 8.314472 / 8.314510
    common_type = promote_type(typeof(ρ), typeof(T), eltype(x))
    _0 = zero(common_type)
    res = _0
    x0 = _0  #for comparison
    ao2 = _0
    ao_zero = _0

    ρc = model.criticalDensity
    Tc = model.criticalTemperature
    nr = model.nr
    zeta = model.zeta
    for i in model.ideal_iters[1]
        x[i] != x0 && begin
            ao2 = ao_zero
            ao3 = ao_zero
            δ = ρ / ρc[i]
            τ = Tc[i] / T
            ao1 = nr[1, i] + nr[2, i] * τ + nr[3, i] * log(τ)
            ao2 =
            nr[4, i] * log(abs(sinh(zeta[1, i] * τ))) -
            nr[5, i] * log(cosh(zeta[2, i] * τ)) +
            nr[6, i] * log(abs(sinh(zeta[3, i] * τ))) -
            nr[7, i] * log(cosh(zeta[4, i] * τ))

            ao3 = log(δ)
            a0 = RR * (ao1 + ao2) + ao3
            res += x[i] * (a0 + log(x[i]))
        end
    end

    for i in model.ideal_iters[2]
        x[i] != x0 && begin
            ao2 = ao_zero
            ao3 = ao_zero
            δ = ρ / ρc[i]
            τ = Tc[i] / T
            ao1 = nr[1, i] + nr[2, i] * τ + nr[3, i] * log(tau)
            ao2 =
                nr[4, i] * log(abs(sinh(zeta[1, i] * tau))) -
                nr[5, i] * log(cosh(zeta[2, i] * tau)) +
                nr[6, i] * log(abs(sinh(zeta[3, i] * tau)))
            ao3 = log(δ)
            a0 = RR * (ao1 + ao2) + ao3
            res += x[i] * (a0 + log(x[i]))
        end
    end

    for i in model.ideal_iters[3]
        x[i] != x0 && begin
            ao2 = ao_zero
            ao3 = ao_zero
            δ = ρ / ρc[i]
            τ = Tc[i] / T
            ao1 = nr[1, i] + nr[2, i] * τ + nr[3, i] * log(τ)
            ao2 =
                nr[4, i] * log(abs(sinh(zeta[1, i] * τ))) -
                nr[5, i] * log(cosh(zeta[2, i] * τ))
            ao3 = log(δ)
            a0 = RR * (ao1 + ao2) + ao3
            res += x[i] * (a0 + log(x[i]))
        end
    end

    for i in model.ideal_iters[4]
        x[i] != x0 && begin
            ao2 = ao_zero
            ao3 = ao_zero
            δ = ρ / ρc[i]
            τ = Tc[i] / T
            ao1 = nr[1, i] + nr[2, i] * τ + nr[3, i] * log(τ)
            ao3 = log(δ)
            a0 = RR * (ao1) + ao3
            res += x[i] * (a0 + log(x[i]))
        end
    end

    return res
end

function __f0(model::GERG2008, rho, T)


    RR = 8.314472 / 8.314510
    common_type = promote_type(typeof(rho), typeof(T))

    res = zero(common_type)
    ao1 = zero(common_type)
    ao2 = zero(common_type)
    ao_zero = zero(common_type)

    delta = rho / model.criticalDensity[1]
    tau = model.criticalTemperature[1] / T

    if !iszero(length(model.ideal_iters[1]))
        ao1 = model.nr[1, 1] + model.nr[2, 1] * tau + model.nr[3, 1] * log(tau)
        ao2 =
            model.nr[4, 1] * log(abs(sinh(model.zeta[1, 1] * tau))) -
            model.nr[5, 1] * log(cosh(model.zeta[2, 1] * tau)) +
            model.nr[6, 1] * log(abs(sinh(model.zeta[3, 1] * tau))) -
            model.nr[7, 1] * log(cosh(model.zeta[4, 1] * tau))

        ao3 = log(delta)
        a0 = RR * (ao1 + ao2) + ao3
        res +=  a0
    end

    if !iszero(length(model.ideal_iters[2]))
        ao1 = model.nr[1, 1] + model.nr[2, 1] * tau + model.nr[3, 1] * log(tau)
        ao2 =
            model.nr[4, 1] * log(abs(sinh(model.zeta[1, 1] * tau))) -
            model.nr[5, 1] * log(cosh(model.zeta[2, 1] * tau)) +
            model.nr[6, 1] * log(abs(sinh(model.zeta[3, 1] * tau)))
        ao3 = log(delta)
        a0 = RR * (ao1 + ao2) + ao3
        res +=  a0
    end

    if !iszero(length(model.ideal_iters[3]))
        ao1 = model.nr[1, 1] + model.nr[2, 1] * tau + model.nr[3, 1] * log(tau)
        ao2 =
            model.nr[4, 1] * log(abs(sinh(model.zeta[1, 1] * tau))) -
            model.nr[5, 1] * log(cosh(model.zeta[2, 1] * tau))
        ao3 = log(delta)
        a0 = RR * (ao1 + ao2) + ao3
        res +=  a0
    end

    if !iszero(length(model.ideal_iters[4]))
        ao1 = model.nr[1, 1] + model.nr[2, 1] * tau + model.nr[3, 1] * log(tau)
        ao3 = log(delta)
        a0 = RR * (ao1) + ao3
        res +=  a0 
    end

    return res
end
