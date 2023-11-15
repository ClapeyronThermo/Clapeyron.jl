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
