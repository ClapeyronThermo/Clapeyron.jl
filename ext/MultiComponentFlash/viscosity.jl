function M.lbc_viscosity(eos::C.EoSModel, p, temperature, ph::M.FlashedPhase{T}; coeff = (0.1023, 0.023364, 0.058533, -0.040758, 0.0093324), shift = -1e-4) where T
    z = ph.mole_fractions

    mw_mix, P_pc, T_pc, V_pc = M.pseudo_critical_properties(eos, z)
    mu_atm = M.atmospheric_mu_estimate(eos, z, temperature)
    e_mix = M.mixture_viscosity_parameter(mw_mix, T_pc, P_pc)
    # From Jossi et al
    # coeffs = [0.1023, 0.023364, 0.058533, -0.040758, 0.0093724, 1e-4]
    # From LBC paper
    # coeffs = [0.1023, 0.023364, 0.058533, -0.040758, 0.0093324, -1e-4]
    # Compute reduced density
    V = M.molar_volume(eos, p, temperature, ph)
    rho_r = V_pc/V
    corr = zero(T)
    for (i, c) in enumerate(coeff)
        corr += c*rho_r^(i-1)
    end
    mu_correction = (corr^4 + shift)/e_mix
    # Final expression is compound and given in centi poise. We convert to Pa s
    # instead for strict SI outputs.
    mu = 1e-3*(mu_atm + mu_correction)
    return mu::T
end

function M.pseudo_critical_properties(model::C.EoSModel, z::V) where V<:AbstractVector{T} where T
    mw = C.molecular_weight(model,z)
    pure = C.split_model(model)
    crit = C.crit_pure.(pure)
    tc = [crit[i][1] for i in 1:length(model)]
    pc = [crit[i][2] for i in 1:length(model)]
    vc = [crit[i][3] for i in 1:length(model)]
    Tc = C.dot(tc,z)
    Pc = C.dot(pc,z)
    Vc = C.dot(vc,z)
    return mw,Pc,Tc,Vc
end

function M.pseudo_critical_properties(model::C.MultiFluid, z::V) where V<:AbstractVector{T} where T
    mw = C.molecular_weight(model,z)
    tc = model.params.Tc.values
    pc = model.params.Pc.values
    vc = model.params.Vc.values
    Tc = C.dot(tc,z)
    Pc = C.dot(pc,z)
    Vc = C.dot(vc,z)
    return mw,Pc,Tc,Vc
end

function M.pseudo_critical_properties(model::C.CubicModel, z::V) where V<:AbstractVector{T} where T
    mw = C.molecular_weight(model,z)
    tc = model.params.Tc.values
    pc = model.params.Pc.values
    Tc = C.dot(tc,z)
    Pc = C.dot(pc,z)
    if hasfield(typeof(model.params),:Vc)
        Vc = C.dot(model.params.Vc.values,z)
    else
        Vc = C.volume(model,Pc,Tc,z)
    end
    return mw,Pc,Tc,Vc
end

function M.atmospheric_mu_estimate(model::C.EoSModel, z::V, temperature) where V<:AbstractVector{T} where T
    a = zero(T)
    b = zero(T)
    Tc = model.params.Tc.values
    Pc = model.params.Pc.values
    Mw = C.mw(model)
    @inbounds for i in 1:length(model)
        mw = Mw[i]/1000.
        T_c = Tc[i]
        p_c = Pc[i]
        zi = z[i]
        # Add contributions to atmospheric mu
        T_r = temperature/T_c
        e_i = M.mixture_viscosity_parameter(mw, T_c, p_c)
        if T_r > 1.5
            mu_i = 17.78e-5*(4.58*T_r - 1.67)^0.625
        else
            mu_i = 34e-5*T_r^(0.94)
        end
        tmp = sqrt(1000.0*mw)*zi
        a += tmp*mu_i/e_i
        b += tmp
    end
    return a/b
end
