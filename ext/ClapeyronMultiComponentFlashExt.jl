module ClapeyronMultiComponentFlashExt
    using Clapeyron: Clapeyron
    using MultiComponentFlash: MultiComponentFlash
    const C = Clapeyron
    const M = MultiComponentFlash
    const S = C.StaticArrays
    const ForwardDiff = C.ForwardDiff

    struct CubicWrapper{T} <: C.EoSModel
        model::T
    end

    #passing GenericCubicEoS to Clapeyron.jl userlocations
    C.can_nt(x::M.GenericCubicEOS) = true
    C.can_nt(x::M.MultiComponentMixture) = true

    function C.to_nt(x::M.MultiComponentMixture)
        res = Dict{Symbol,Any}()
        if x.binary_interaction !== nothing
            res[:k] = x.binary_interaction
        end
        props = x.properties
        res[:Mw] = [1000*getfield(props_i,:mw) for props_i in props]
        res[:Tc] = [getfield(props_i,:T_c) for props_i in props]
        res[:Pc] = [getfield(props_i,:p_c) for props_i in props]
        res[:Vc] = [getfield(props_i,:V_c) for props_i in props]
        res[:acetricfactor] = [getfield(props_i,:ω) for props_i in props]
        return res
    end
    function C.to_nt(x::M.GenericCubicEOS)
        res = C.to_nt(x.mixture)
        #TODO:check that volume shift
    end
    M.number_of_components(model::C.EoSModel) = length(model)

    function M.mixture_compressibility_factor(eos::C.EoSModel, cond,
        forces = force_coefficients(eos, cond),
        scalars = force_scalars(eos, cond, forces))
        C.compressibility_factor(eos,cond.p,cond.T,cond.z)
    end

    function M.molar_volume(model::C.EoSModel, p, T, ph::M.FlashedPhase{X}) where X
        convert(X,C.volume(model,p,T,ph.Z))
    end

    function M.mass_density(model::C.EoSModel, p, T, ph::M.FlashedPhase{X}) where X
        convert(X,C.mass_density(model,p,T,ph.Z))
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

    function M.initial_guess_K!(K, eos::C.EoSModel, cond)
        C.tp_flash_K0!(K, eos, cond.p, cond.T)
    end

    function M.flash_storage_internal!(out, eos::C.EoSModel, cond, method; inc_jac = isa(method, M.AbstractNewtonFlash), static_size = false, kwarg...)
        n = M.number_of_components(eos)
        TT = typeof(cond.p+cond.T+one(eltype(eos)))
        splt = C.split_model(eos)
        out[:split_model] = splt
        out[:forces] = nothing
        out[:crit] = C.crit_pure.(splt)
        if static_size
            alloc_vec = () -> C.StaticArrays.@MVector zeros(TT,n)
        else
            alloc_vec = () -> zeros(TT,n)
        end
        out[:x] = alloc_vec()
        out[:y] = alloc_vec()
        #for fugacities
        out[:buffer1] = alloc_vec()
        out[:buffer2] = alloc_vec()
        if inc_jac
            M.flash_storage_internal_newton!(out, eos, cond, method, static_size = static_size; kwarg...)
        end
        return out
    end

    #this is only defined with cubic EoS.
    M.force_coefficients(eos::C.EoSModel, cond) = nothing
    M.force_scalars(eos::C.EoSModel, cond, forces) = nothing
    M.force_coefficients!(forces, eos::C.EoSModel, c) = nothing

    function M.mixture_fugacities!(f, eos::C.EoSModel, cond, forces = M.force_coefficients(eos, cond), scalars = M.force_scalars(eos, cond, forces))
        phase::Symbol = get(cond,:phase,:unknown)
        C.fugacity_coefficient!(f,eos,cond.p,cond.T,cond.z,phase = phase)
        f .*= cond.p
        f .= f .* cond.z
        return f
    end

    #TODO: this could be removed
    function M.stability_2ph!(storage, K, eos::C.EoSModel, c; verbose = false, kwarg...)
        forces = storage.forces
        f_z = storage.buffer1
        f_xy = storage.buffer2
        x, y = storage.x, storage.y
        z, p, T = c.z, c.p, c.T
        liquid = (p = p, T = T, z = x, phase = :liquid)
        vapor = (p = p, T = T, z = y,phase = :vapor)
        # Update fugacities for current conditions used in both tests
        M.mixture_fugacities!(f_z, eos, c, forces)
        #props = eos.mixture.properties
        C.wilson_k_values!(K,eos,p,T,storage.crit)
        stable_vapor, i_v = M.michelsen_test!(vapor, f_z, f_xy, vapor.z, z, K, eos, c, forces, Val(true); kwarg...)
        C.wilson_k_values!(K,eos,p,T,storage.crit)
        stable_liquid, i_l = M.michelsen_test!(liquid, f_z, f_xy, liquid.z, z, K, eos, c, forces, Val(false); kwarg...)
        stable = stable_vapor && stable_liquid
        if !stable
            @. K = y/x
        end
        if verbose
            @info "Stability done. Iterations:\nV: $i_v\nL: $i_l" stable_vapor stable_liquid stable
        end
        return stable
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

    function M.flash_update!(K, storage, type::M.SSIFlash, eos::C.EoSModel, cond, forces, V::F, iteration) where F
        z = cond.z
        x, y = storage.x, storage.y
        p, T = cond.p, cond.T
        x = M.liquid_mole_fraction!(x, z, K, V)
        y = M.vapor_mole_fraction!(y, x, K)
        phi_l = storage.buffer1
        phi_v = storage.buffer2
        liquid = (p = p, T = T, z = x,phase = :liquid)
        vapor = (p = p, T = T, z = y,phase = :vapor)
        ϵ = zero(F)
        M.mixture_fugacities!(phi_l, eos, liquid, forces)
        M.mixture_fugacities!(phi_v, eos, vapor, forces)
        r = phi_l
        r ./= phi_v
        K .*= r
        ϵ = mapreduce(ri -> abs(1-ri), max ,r)
        V = C.rachfordrice(K, z; β0=V)
        return (V, ϵ)::Tuple{F, F}
    end

    #support for tp_flash interface

    function C.MCFlashJL(;method = M.SSIFlash(),
        storage = nothing,
        V = NaN,
        K = nothing,
        tolerance = 1e-8,
        maxiter = 20000,
        verbose = false,
        check = true,
        z_min = nothing,
        update_forces = true)
        kwargs = (;tolerance,maxiter,verbose,check,z_min,update_forces)
        return C.MCFlashJL(method,storage,K,V,kwargs)
    end

    function M.flash_storage(model::C.EoSModel,p,T,z,method::C.MCFlashJL)
        cond = (p = p,T = T,z = z)
        mcf_method = method.method
        M.flash_storage(model,cond;method = mcf_method,method.kwarg...)
    end

    function C.tp_flash_impl(model,p,T,z,method::C.MCFlashJL)
        S = method.storage == nothing ? M.flash_storage(model,p,T,z,method) : method.storage
        K = method.K === nothing ? C.wilson_k_values!(zeros(typeof(p+T+one(eltype(model))),length(model)),model,p,T,S.crit) : method.K

        conditions = (p = p, T = T,z = z)
        V = M.flash_2ph!(S,K,model,conditions;method = method.method,method.kwarg...)
        x = M.liquid_mole_fraction!(S.x, z, K, V)
        y = M.vapor_mole_fraction!(S.y, x, K)
        G = C.__tpflash_gibbs_reduced(model,p,T,x,y,V,:vle)
        X = hcat(x,y)'
        nvals = X.*[1-V
                V] .* sum(z)
        return (X, nvals, G)
    end

    #flow_coupler.jl
    function M.flashed_mixture_2ph(eos::C.EoSModel, cond, _K = nothing; method = M.SSIFlash(), kwarg...)
        # Convenience function for getting flashed phases
        S = M.flash_storage(eos, cond, method = method)
        p,T,z = cond.p,cond.T,cond.z
        K = _K === nothing ? C.wilson_k_values!(zeros(typeof(p+T+one(eltype(model))),length(model)),model,p,T,S.crit) : _K
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
            state = two_phase_lv
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
end #module

