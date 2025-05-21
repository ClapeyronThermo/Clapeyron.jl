function M.flash_storage_internal!(out, eos::C.EoSModel, cond, method; inc_jac = isa(method, M.AbstractNewtonFlash), static_size = false,inc_bypass = false, kwarg...)
    n = M.number_of_components(eos)
    TT = typeof(one(eltype(eos)))
    splt = C.split_pure_model(eos)
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
    out[:buffer1] = alloc_vec()
    out[:buffer2] = alloc_vec()
    if inc_jac
        M.flash_storage_internal_newton!(out, eos, cond, method, static_size = static_size; kwarg...)
    end

    if inc_bypass
        out[:bypass] = michelsen_critical_point_measure_storage(eos, static_size = static_size)
    end

    return out
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

function M.update_flash_jacobian!(J, r, eos::C.EoSModel, p, T, z, x, y, V, forces)
    has_r = !isnothing(r)
    n = M.number_of_components(eos)
    liquid = (p = p, T = T, z = x)
    vapor = (p = p, T = T, z = y)

    #those 3 lines are the only difference. we need an overload here
    Z_l = C.compressibility_factor(eos,p,T,x,phase = :l)
    Z_v = C.compressibility_factor(eos,p,T,y,phase = :v)
    s_l,s_v = nothing,nothing

    if isa(V, ForwardDiff.Dual)
        np = length(V.partials)
    else
        np = length(p.partials)
    end
    # Isofugacity constraint
    # f_li - f_vi ∀ i
    @inbounds for c in 1:n
        f_l = M.component_fugacity(eos, liquid, c, Z_l, forces, s_l)
        f_v = M.component_fugacity(eos, vapor, c, Z_v, forces, s_v)
        Δf = f_l - f_v
        if has_r
            @inbounds r[c+n] = Δf.value
        end
        for i = 1:np
            @inbounds J[c+n, i] = Δf.partials[i]
        end
    end
    # x_i*(1-V) - V*y_i - z_i = 0 ∀ i
    # Σ_i x_i - y_i = 0
    T = Base.promote_type(eltype(x), eltype(y))
    L = 1 - V
    Σxy = zero(T)
    @inbounds for c in 1:n
        xc, yc = x[c], y[c]
        Σxy += (xc - yc)
        M = L*xc + V*yc - z[c]
        if has_r
            @inbounds r[c] = M.value
        end
        @inbounds for i = 1:np
            J[c, i] = MultiComponentFlash.∂(M, i)
        end
    end
    if has_r
        @inbounds r[end] = Σxy.value
    end
    if isa(Σxy, ForwardDiff.Dual)
        @inbounds for i = 1:np
            J[end, i] = M.∂(Σxy, i)
        end
    else
        @. J[end, :] = 0
    end
end
