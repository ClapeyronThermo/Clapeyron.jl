function M.FlashedMixture2Phase(eos::C.EoSModel, T = Float64, T_num = Float64)
    n = M.number_of_components(eos)
    V = zero(T)
    # K values are always doubles
    K = zeros(T_num, n)
    liquid = M.FlashedPhase(n, T)
    vapor = M.FlashedPhase(n, T)
    return M.FlashedMixture2Phase(M.unknown_phase_state_lv, K, V, liquid, vapor)
end

function M.single_phase_label(model::C.EoSModel,cond)
    V,_,_ = C._label_and_volumes(model,cond)
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
    x = similar(K)
    y = similar(K)
    if V == 0 || V == 1 || isnan(V)
        #=
        instead of Li correlation, we just check the gibbs energies of the mixtures.
        we do this anyway in the volume calculation, so it better to be more explicit.
        =#  
        V,Vl,Vv = C._label_and_volumes(eos, conditions)
        n = sum(z)
        nRT = C.Rgas(eos)*T*n
        if V == 0
            state = M.single_phase_l
            Zx = Vl*p/nRT
        else
            state = M.single_phase_v
            Zx = Vv*p/nRT
        end
        Z_L,Z_V = Zx,Zx
        x .= z ./ n
        y .= z ./ n
    else
        @. x = M.liquid_mole_fraction(z, K, V)
        @. y = M.vapor_mole_fraction(x, K)
        state = M.two_phase_lv
        Z_L = C.compressibility_factor(eos, p, T, x, phase = :l)
        Z_V = C.compressibility_factor(eos, p, T, y, phase = :v)
    end
    nan = zero(V)/zero(V)
    return M.FlashedMixture2Phase(state, K, V, x, y, Z_L, Z_V, nan, conditions, rep.stability)
end

@inline function M.two_phase_vapor_saturation(eos::C.EoSModel, p, Temp, f::M.FlashedMixture2Phase{T}) where T
    state = f.state
    # A faster definition that doesn't go via molar volume, but assumes no shifts:
    V = f.V
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