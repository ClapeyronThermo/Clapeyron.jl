function M.FlashedMixture2Phase(eos::C.EoSModel, T = Float64, T_num = Float64)
    n = M.number_of_components(eos)
    V = zero(T)
    # K values are always doubles
    K = zeros(T_num, n)
    liquid = M.FlashedPhase(n, T)
    vapor = M.FlashedPhase(n, T)

    return M.FlashedMixture2Phase(M.unknown_phase_state_lv, K, V, liquid, vapor)
end

function _label_and_volumes(model::C.EoSModel,cond)
    #gibbs comparison, the phase with the least amount of gibbs energy is the most stable.
    p,T,z = cond.p,cond.T,cond.z
    Vl = C.volume(model,p,T,z,phase =:l)
    Vv = C.volume(model,p,T,z,phase =:v)
    function gibbs(fV)
        isnan(fV) && return one(fV)/zero(fV)
        _df,_f =  C.∂f(model,fV,T,z)
        dV,_ = _df
        return ifelse(abs((p+dV)/p) > 0.03,zero(dV)/one(dV),_f + p*fV)
    end
    isnan(Vl) && return 1,Vv,Vv #could not converge on gas volume, assuming stable liquid phase
    isnan(Vv) && return 0,Vl,Vl #could not converge on liquid volume, assuming stable gas phase
    gl,gv = gibbs(Vl),gibbs(Vv)
    V = gv < gl ? 1 : 0
    return V,Vl,Vv
end

function M.single_phase_label(model::C.EoSModel,cond)
    V,_,_ = _label_and_volumes(model,cond)
    return V
end



function M.flashed_mixture_2ph(eos::C.EoSModel, cond, _K = nothing; method = M.SSIFlash(), kwarg...)
    # Convenience function for getting flashed phases
    S = M.flash_storage(eos, cond, method = method)
    p,T,z = cond.p,cond.T,cond.z
    K = _K === nothing ? C.wilson_k_values!(zeros(typeof(p+T+one(eltype(eos))),length(eos)),eos,p,T,S.crit) : _K
    return M.flashed_mixture_2ph!(S, eos, cond, K; method = method, kwarg...)
end

function M.flashed_mixture_2ph!(storage, eos::C.EoSModel, conditions, K; kwarg...)
    V, K, rep = M.flash_2ph!(storage, K, eos, conditions; kwarg..., extra_out = true)
    p = conditions.p
    T = conditions.T
    z = conditions.z

    x = @. M.liquid_mole_fraction(z, K, V)
    y = @. M.vapor_mole_fraction(x, K)
    if V == 0 || V == 1
    #=
    instead of Li correlation, we just check the gibbs energies of the mixtures.
    we do this anyway in the volume calculation, so it better to be more explicit.
    =#
        V,Vl,Vv = _label_and_volumes(eos, conditions)
        if V == 0
            state = single_phase_l
            Zx = Vl*p/C.Rgas(model)/T/sum(z)
        else
            state = single_phase_v
            Zx = Vl*p/C.Rgas(model)/T/sum(z)
        end
        Z_L,Z_V = Zx,Zx
    else
        state = M.two_phase_lv
        Z_L = C.compressibility_factor(eos, p, T, x, phase = :l)
        Z_V = C.compressibility_factor(eos, p, T, y, phase = :v)
    end
    return M.FlashedMixture2Phase(state, K, V, x, y, Z_L, Z_V)
end

@inline function M.two_phase_vapor_saturation(eos::C.EoSModel, p, Temp, f::M.FlashedMixture2Phase{T}) where T
    state = f.state
    # A faster definition that doesn't go via molar volume, but assumes no shifts:
    Z_l = f.liquid.Z
    Z_v = f.vapor.Z
    S_v = Z_v*V/(Z_l*(1-V) + Z_v*V)
    return S_v
    #the volume calculated in
    #=
    V = f.V
    L = one(V) - V
    vol_v = V*molar_volume(eos, p, Temp, f.vapor)
    vol_l = L*molar_volume(eos, p, Temp, f.liquid)
    S_v = vol_v/(vol_v + vol_l)
    return S_v =#
end

function M.molar_volume(model::C.EoSModel, p, T, ph::M.FlashedPhase{X}) where X
    Z = ph.Z
    v = convert(X,Z*C.Rgas(model)*T/p)
    return v
end

function M.mass_density(model::C.EoSModel, p, T, ph::M.FlashedPhase{X}) where X
    Z = ph.Z
    v = Z*C.Rgas(model)*T/p
    molar_weight = C.molecular_weight(model,ph.mole_fractions)
    return convert(X,molar_weight/v)
end

function mass_densities(model::C.EoSModel, p, T, f::M.FlashedMixture2Phase{X}) where X
    state = f.state
    # @assert state != unknown_phase_state_lv "Phase state is not known. Cannot compute densities. Has flash been called?."
    if M.liquid_phase_present(state)
        l = C.mass_density(model,p,T,f.liquid,phase =:l)::X
    else
        l = zero(X)
    end
    if M.vapor_phase_present(state)
        v = C.mass_density(model,p,T,f.vapor,phase = :v)::X
    else
        v = zero(X)
    end
    return (ρ_l = l::X, ρ_v = v::X)
end
