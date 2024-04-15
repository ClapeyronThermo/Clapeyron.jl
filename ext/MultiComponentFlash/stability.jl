    #TODO: this could be removed
    function M.stability_2ph!(storage, K, eos::C.EoSModel, c;
        verbose::Bool = false,
        extra_out::Bool = false,
        check_vapor::Bool = true,
        check_liquid::Bool = true,
        kwarg...)
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
        
        if check_vapor
            C.wilson_k_values!(K,eos,p,T,storage.crit)
            vv = M.michelsen_test!(vapor, f_z, f_xy, vapor.z, z, K, eos, c, forces, Val(true); kwarg...)
        else
            vv = (true,true,0)
        end

        if check_vapor
            C.wilson_k_values!(K,eos,p,T,storage.crit)
            ll = M.michelsen_test!(liquid, f_z, f_xy, liquid.z, z, K, eos, c, forces, Val(false); kwarg...)
        else
            ll = (true,true,0)
        end
        
        stable_vapor, trivial_vapor, i_v = vv
        stable_liquid, trivial_liquid, i_l = ll

        report = M.StabilityReport(
            stable_liquid = stable_liquid,
            trivial_liquid = trivial_vapor,
            stable_vapor = stable_vapor,
            trivial_vapor = trivial_vapor
        )

        stable = report.stable
        if !stable
            @. K = y/x
        end
        if verbose
            @info "Stability done. Iterations:\nV: $i_v\nL: $i_l" stable_vapor stable_liquid stable
        end
        if extra_out
            out = (stable, report)
        else
            out = stable
        end
        return out
    end
