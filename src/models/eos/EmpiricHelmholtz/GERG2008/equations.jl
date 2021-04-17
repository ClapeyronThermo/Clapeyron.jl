


_gerg_asymetric_mix_rule(xi, xj, b) = b * (xi + xj) / (xi * b^2 + xj)

function mol_volume(mixmodel::MixRule{GERG2008{MULTI},X},::CriticalPoint,unit=u"m^3/mol") where X
    x = mol_fraction(mixmodel) 
    model = thermomodel(mixmodel)
    vc = mixing_rule_asymetric(
        (a,b) -> ((cbrt(a) + cbrt(b))*0.5)^3,
        _gerg_asymetric_mix_rule,
        x,
        model.criticalVolume,
        model.gamma_v,
        model.beta_v,
    )
    return convert_unit(u"m^3/mol",unit,vc)
end 

function _delta(model::GERG2008{MULTI}, rho, T, x)
    vcmix = mol_volume(MixRule(model,x),CriticalPoint())
    return rho * vcmix
end

#==
function _tau(model::GERG2008{MULTI}, rho, T, x)
    Tr = mixing_rule_asymetric(
        geometric_mean_rule,
        _gerg_asymetric_mix_rule,
        x,
        model.criticalTemperature,
        model.gamma_T,
        model.beta_T,
    )
    return Tr / T
end
==#

function _tau(model::GERG2008{MULTI}, rho, T, x)
    Tcmix = temperature(MixRule(model,x),CriticalPoint())
    return Tcmix / T
end

function temperature(mixmodel::MixRule{GERG2008{MULTI},X},::CriticalPoint,unit=u"K") where {X}
    x = mol_fraction(mixmodel) 
    model = thermomodel(mixmodel)
    Tcmix = mixing_rule_asymetric(
        (a,b)->sqrt(a * b),
        _gerg_asymetric_mix_rule,
        x,
        model.criticalTemperature,
        model.gamma_T,
        model.beta_T,
    )
    return convert_unit(u"K",unit,Tcmix)
end 


function _delta(model::GERG2008{SINGLE}, rho, T)
    rhor = inv(only(model.criticalVolume))
    return rho / rhor
end

function _tau(model::GERG2008{SINGLE}, rho, T)
    Tr = only(model.criticalTemperature)
    return Tr / T
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

    for i = 1:length(model.syms)
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



@inline function mol_helmholtz0_impl(mt::MultiVT,model::GERG2008{MULTI}, v, T, x)
    rho = 1.0e-3 / v
    R = RGAS
    return R * T * _f0(model, rho, T, x)
end

@inline function mol_helmholtz0_impl(mt::SingleVT,model::GERG2008{SINGLE}, v, T)
    rho = 1.0e-3 / v
    R = RGAS
    return R * T * _f0(model, rho, T)
end

@inline function αR_impl(mt::MultiVT,model::GERG2008{MULTI}, _rho, T, x)
    rho = 1.0e-3 * _rho
    R = RGAS
    delta = _delta(model, rho, T, x)
    tau = _tau(model, rho, T, x)
    return _fr1(model, delta, tau, x) + _fr2(model, delta, tau, x)
end

@inline function αR_impl(mt::SingleVT,model::GERG2008{SINGLE}, _rho, T)
    rho = 1.0e-3 * _rho
    R = RGAS
    delta = _delta(model, rho, T)
    tau = _tau(model, rho, T)
    return _fr1(model, delta, tau)
end


@inline function mol_helmholtzR_impl(mt::MultiVT,model::GERG2008{MULTI}, v, T, x)
    rho = 1.0e-3 / v
    R = RGAS
    delta = _delta(model, rho, T, x)
    tau = _tau(model, rho, T, x)
    return R * T * (_fr1(model, delta, tau, x) + _fr2(model, delta, tau, x))
end

@inline function mol_helmholtzR_impl(mt::SingleVT,model::GERG2008{SINGLE}, v, T)
    rho = 1.0e-3 / v
    R = RGAS
    delta = _delta(model, rho, T)
    tau = _tau(model, rho, T)
    return R * T * _fr1(model, delta, tau)
end


@inline function mol_helmholtz_impl(mt::MultiVT,model::GERG2008{MULTI}, v, T, x)
    rho = 1.0e-3 / v
    R = RGAS
    delta = _delta(model, rho, T, x)
    tau = _tau(model, rho, T, x)
    return R *
           T *
           (
               _f0(model, rho, T, x) +
               _fr1(model, delta, tau, x) +
               _fr2(model, delta, tau, x)
           )
end


@inline function mol_helmholtz_impl(mt::SingleVT,model::GERG2008{SINGLE}, v, T)
    rho = 1.0e-3 / v
    R = RGAS
    delta = _delta(model, rho, T)
    tau = _tau(model, rho, T)
    return R *
           T *
           (
               _f0(model, rho, T) +
               _fr1(model, delta, tau)
           )
end



molecular_weight(model::GERG2008{MULTI}) = model.molecularWeight
molecular_weight(model::GERG2008{SINGLE}) = only(model.molecularWeight)



#critical, single component:
mol_density(model::GERG2008{SINGLE},::CriticalPoint,unit=u"mol/(m^3)") = convert_unit(u"mol/L",unit,only(model.criticalDensity))
pressure(model::GERG2008{SINGLE},::CriticalPoint,unit=u"Pa") = convert_unit(u"Pa",unit,only(model.criticalPressure))
temperature(model::GERG2008{SINGLE},::CriticalPoint,unit=u"K") = convert_unit(u"K",unit,only(model.criticalTemperature))
mol_volume(model::GERG2008{SINGLE},::CriticalPoint,unit=u"m^3/mol") = convert_unit(u"m^3/mol",unit,only(model.criticalVolume))
acentric_factor(model::GERG2008{SINGLE}) = only(model.acentric_factor)

#critical, multi component:

function mol_density(model::GERG2008{MULTI},::CriticalPoint,unit=u"mol/(m^3)")
    return convert_unit.(u"mol/L",unit,model.criticalDensity)
end

function pressure(model::GERG2008{MULTI},::CriticalPoint,unit=u"Pa")
    return convert_unit.(u"Pa",unit,model.criticalPressure)
end

function temperature(model::GERG2008{MULTI},::CriticalPoint,unit=u"K")
    return convert_unit.(u"K",unit,model.criticalTemperature)
end

function mol_volume(model::GERG2008{MULTI},::CriticalPoint,unit=u"m^3/mol")
    return convert_unit.(u"m^3/mol",unit,model.criticalVolume)
end

function acentric_factor(model::GERG2008{MULTI})
    return model.acentric_factor
end


    

single_sat_aprox(model::GERG2008{SINGLE}) = LeeKesler(model)
volume_solver_type(model::GERG2008) = SUVA()

function Base.show(io::IO, sp::GERG2008{MULTI})
    ln = length(sp.syms)
    println(io,ln,"-element GERG008 helmholtz model, with compounds:")
    for i in 1:ln-1
        println(io," n",unicode_subscript(i)," : ",sp.syms[i])
    end
    print(io," n",unicode_subscript(ln)," : ",sp.syms[ln])
end

function Base.show(io::IO, sp::GERG2008{SINGLE})
    compound = only(sp.syms)
    print(io,"GERG008 helmholtz model for ",compound)
end

export GERG2008
 