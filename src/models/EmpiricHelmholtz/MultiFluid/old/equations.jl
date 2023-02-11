
#=

function T_scale(model::GERG2008,x=SA[1.])
    return mixing_rule_asymetric(
        (a,b)->sqrt(a * b),
        _gerg_asymetric_mix_rule,
        x,
        model.Tc,
        model.gamma_T,
        model.beta_T,
    )
end

# function T_scale(model::GERG2008)
#     return only(model.Tc)
# end

function T_scales(model::GERG2008,x=SA[1.])
    return model.Tc
end

# function p_scale(model::GERG2008)
#     return only(model.Pc)
# end

function p_scale(model::GERG2008,x=SA[1.])
    return dot(x,model.Pc)
end

# function p_scales(model::GERG2008,x=SA[1.])
#     return only(model.Pc)
# end

function _v_scale(model::GERG2008,x=SA[1.])
    return mixing_rule_asymetric(
        (a,b) -> ((cbrt(a) + cbrt(b))*0.5)^3,
        _gerg_asymetric_mix_rule,
        x,
        model.vc,
        model.gamma_v,
        model.beta_v,
    )
end

# function _v_scale(model::GERG2008)
#     return only(model.vc)
# end

# function lb_volume(model::GERG2008)
#     return only(model.lb_v)
# end

function lb_volume(model::GERG2008,x=SA[1.];phase=:unknown)
    return dot(x,model.lb_v)
end

function _delta(model::GERG2008, rho, T, x=SA[1.])
    vcmix = _v_scale(model,x)
    return rho * vcmix
end

function _tau(model::GERG2008, rho, T, x=SA[1.])
    Tcmix  = T_scale(model,x)
    return Tcmix / T
end

# function _delta(model::GERG2008, rho, T)
#     vcmix = _v_scale(model)
#     return rho * vcmix
# end

# function _tau(model::GERG2008, rho, T)
#     Tcmix  = T_scale(model)
#     return Tcmix / T
# end

function _f0(model::GERG2008, ρ, T, x)
    RR = 8.314472 / 8.314510
    common_type = promote_type(typeof(ρ), typeof(T), eltype(x))
    _0 = zero(common_type)
    res = _0
    x0 = _0  #for comparison
    ao2 = _0
    ao_zero = _0

    ρc = model.rhoc
    Tc = model.Tc
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
            ao1 = nr[1, i] + nr[2, i] * τ + nr[3, i] * log(τ)
            ao2 =
                nr[4, i] * log(abs(sinh(zeta[1, i] * τ))) -
                nr[5, i] * log(cosh(zeta[2, i] * τ)) +
                nr[6, i] * log(abs(sinh(zeta[3, i] * τ)))
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

function _f0(model::GERG2008, rho, T)


    RR = 8.314472 / 8.314510
    common_type = promote_type(typeof(rho), typeof(T))

    res = zero(common_type)
    ao1 = zero(common_type)
    ao2 = zero(common_type)
    ao_zero = zero(common_type)

    delta = rho / model.rhoc[1]
    tau = model.Tc[1] / T

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


function _fr1(model::GERG2008, delta, tau) 
    common_type = promote_type(typeof(delta), typeof(tau))
    res = zero(common_type)
    res0 = zero(common_type)
    res1 = zero(common_type)

    res1 = res0
    for k = 1:model.k_pol_ik[1]
        res1 += model.n0ik[1][k] * (delta^model.d0ik[1][k]) * (tau^model.t0ik[1][k])
    end

    for k = (model.k_pol_ik[1]+1):(model.k_exp_ik[1]+model.k_pol_ik[1])
        res1 +=
            model.n0ik[1][k] *
            (delta^model.d0ik[1][k]) *
            (tau^model.t0ik[1][k]) *
            exp(-delta^model.c0ik[1][k-model.k_pol_ik[1]])
    end
    res += res1
    return res
end

function _fr1(model::GERG2008, delta, tau,x)
    common_type = promote_type(typeof(delta), typeof(tau), eltype(x))
    res = zero(common_type)
    res0 = zero(common_type)
    res1 = zero(common_type)
    x0 = zero(common_type)

    for i = 1:length(model.components)
        x[i] != x0 && begin
            res1 = res0
            for k = 1:model.k_pol_ik[i]
                res1 += model.n0ik[i][k] * (delta^model.d0ik[i][k]) * (tau^model.t0ik[i][k])
            end

            for k = (model.k_pol_ik[i]+1):(model.k_exp_ik[i]+model.k_pol_ik[i])
                res1 +=
                    model.n0ik[i][k] *
                    (delta^model.d0ik[i][k]) *
                    (tau^model.t0ik[i][k]) *
                    exp(-delta^model.c0ik[i][k-model.k_pol_ik[i]])
            end
            res += x[i] * res1

        end
    end

    return res
end




function _fr2(model::GERG2008, delta, tau, x)
    common_type = promote_type(typeof(delta), typeof(tau), eltype(x))
    res = zero(common_type)
    res0 = zero(common_type)
    x0 = zero(common_type)

    for kk = 1:length(model.Aij_indices)
        i1, i2, i0 = model.Aij_indices[kk]
        x[i1] != x0 && x[i2] != x0 && begin
            res1 = res0
            for j = 1:model.k_pol_ijk[i0]
                res1 +=
                    model.nijk[i0][j] * (delta^model.dijk[i0][j]) * (tau^model.tijk[i0][j])
            end

            for j = (model.k_pol_ijk[i0]+1):(model.k_pol_ijk[i0]+model.k_exp_ijk[i0])


                idx = j - model.k_pol_ijk[i0]

                res1 +=
                    model.nijk[i0][j] *
                    (delta^model.dijk[i0][j]) *
                    (tau^model.tijk[i0][j]) *
                    exp(
                        -model.etaijk[i0][idx] * (delta - model.epsijk[i0][idx])^2 -
                        model.betaijk[i0][idx] * (delta - model.gammaijk[i0][idx]),
                    )
            end
            res += res1 * x[i1] * x[i2] * model.fij[i0]

        end
    end
    return res
end

function a_ideal(model::GERG2008, v, T, z=SA[1.0])
    N = sum(z)
    len = length(z)
    x = z/N
    rho = 1.0e-3 / (v/N)
    if len == 1
        return _f0(model, rho, T) 
    else
        return _f0(model, rho, T, x)
    end
end

function a_res(model::GERG2008, v, T, z=SA[1.0])
    N = sum(z)
    len = length(z)
    rho = 1.0e-3 / (v/N)
    x = z/N
    delta = _delta(model, rho, T, x)
    tau = _tau(model, rho, T, x)
    if len == 1
        return _fr1(model, delta, tau) 
    else
        return (_fr1(model, delta, tau, x) + _fr2(model, delta, tau, x))
    end
end

function eos(model::GERG2008, v, T, z=SA[1.0])
    negative_vt(v,T) && return nan_num(v,T,z)
    N = sum(z)
    len = length(z)
    rho = 1.0e-3 / (v/N)
    x = z/N
    R = R̄
    delta = _delta(model, rho, T, x)
    tau = _tau(model, rho, T, x)
    if len == 1
        return N * R * T * (_f0(model, rho, T) + _fr1(model, delta, tau)) 
    else
        return N * R * T * (_f0(model, rho, T, x)+ _fr1(model, delta, tau, x) + _fr2(model, delta, tau, x))
    end
end
=#
