# Fixed-composition criticality: concentration scaling and a nonzero cofactor pivot.
# Assumes a simple zero eigenvalue; index_reduction has removed zero amounts.
function scaled_mixture_critical_constraint(model,V,T,z)
    fr(x) = sum(x)*a_res(model,V,T,x)
    function B(x)
        sx = sqrt.(x)
        BB = ForwardDiff.hessian(fr,x) .* (sx * transpose(sx))
        for i in eachindex(x)
            BB[i,i] += one(eltype(BB))
        end
        return BB
    end
    L(x) = det(Symmetric(B(x)))
    HH = B(z)
    LL = det(Symmetric(HH))
    dL = ForwardDiff.gradient(L,z) .* sqrt.(z)
    # At a simple critical point, principal cofactors are proportional to q_i^2.
    # Select a nonvanishing null-vector coordinate rather than always the last one.
    ids = collect(eachindex(z))
    cofactors = map(ids) do i
        keep = filter(!=(i),ids)
        det(Symmetric(HH[keep,keep]))
    end
    pivot = argmax(abs.(cofactors))
    HH[pivot,:] .= dL
    return LL,det(HH)
end

function Obj_crit_mix(model::EoSModel,F,z,V,T)
    L,detM = scaled_mixture_critical_constraint(model,V,T,z)
    F[1] = L
    F[2] = detM
    return F
end


"""
    crit_mix(model::EoSModel,z;v0=x=x0_crit_mix(model,z))

Returns the critical mixture point at a given composition.

Returns a tuple, containing:
- Critical Mixture Temperature `[K]`
- Critical Mixture Pressure `[Pa]`
- Critical Mixture Volume `[m³]`
"""
function crit_mix(model::EoSModel,z;v0=nothing)
    ∑z = sum(z)
    model_r,idx_r = index_reduction(model,z)

    if length(model_r)==1
        (T_c,p_c,V_c) = crit_pure(model_r)
        return (T_c,p_c,V_c*∑z)
    end


    z_r = z[idx_r]
    z_r ./= ∑z
    if v0 === nothing
        v0 = x0_crit_mix(model_r,z_r)
    end

    x0 = [v0[1],v0[2]] #could replace for MVector{2}
    f! = (F,x) -> Obj_crit_mix(model_r, F, z_r, exp10(x[1]), x[2])
    r  = Solvers.nlsolve(f!,x0,LineSearch(Newton2(x0),Backtracking()))
    sol = Solvers.x_sol(r)
    !__check_convergence(r) && (sol .= NaN)
    T_c = sol[2]
    V_c = exp10(sol[1])
    p_c = pressure(model_r, V_c, T_c, z_r)
    return (T_c, p_c, ∑z*V_c)
end


"""
    x0_crit_mix(model::EoSModel,z)

Initial point for `crit_mix(model,z)`.

Returns a tuple, containing:
- Base 10 logarithm of initial guess for critical molar Volume `[m³·mol⁻¹]`
- Initial guess for critical temperature `[K]`

If the model is not splittable or is it an electrolyte model (where the ions or salts cannot be considered as standalone components), then a mechanical critical point approximation is calculated first via [`mechanical_crit_bracket`](@ref) 
"""
function x0_crit_mix(model::EoSModel,z)
    if is_splittable(model)
        pure = split_pure_model(model)
        crit = crit_pure.(pure)
        vci = getindex.(crit,3)
        tci = getindex.(crit,1)
        ∑z = sum(z)
        V_c = dot(z,vci)/∑z
        T_c  = prod(tci[i]^(z[i]/∑z) for i ∈ 1:length(model))
        return (log10(V_c),T_c)
    else
        return x0_crit_bisection(model,z)
    end
end
