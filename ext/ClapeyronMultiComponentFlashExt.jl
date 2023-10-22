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

    C.can_nt(x::M.GenericCubicEOS) = true
    C.to_nt(x::M.GenericCubicEOS) = C.to_nt(x.mixture)
    
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

    M.number_of_components(model::C.EoSModel) = length(model)

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
        out[:phi_l] = alloc_vec()
        out[:phi_v] = alloc_vec()

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


    function M.flash_update!(K, storage, type::M.SSIFlash, eos::C.EoSModel, cond, forces, V::F, iteration) where F
        z = cond.z
        x, y = storage.x, storage.y
        p, T = cond.p, cond.T
        x = M.liquid_mole_fraction!(x, z, K, V)
        y = M.vapor_mole_fraction!(y, x, K)
        phi_l = storage.phi_l
        phi_v = storage.phi_v
        liquid = (p = p, T = T, z = x,phase = :liquid)
        vapor = (p = p, T = T, z = y,phase = :vapor)
        ϵ = zero(F)
        M.mixture_fugacities!(phi_l, eos, liquid, forces)
        M.mixture_fugacities!(phi_v, eos, vapor, forces)
        r = phi_l
        r ./= phi_v
        K .*= r
        ϵ = mapreduce(ri -> abs(1-ri), max ,r)
        return (V, ϵ)::Tuple{F, F}
    end
end #module