"""
    crit_pure(model::EoSModel,x0=nothing)

Calculates the critical point of a single component modelled by `model`. 
Returns a tuple, containing:

- Critical temperature `[K]`
- Critical pressure `[Pa]`
- Critical molar volume `[m³·mol⁻¹]`
"""
crit_pure

function crit_pure(model::EoSModel)
    satmodel = saturation_model(model)
    if satmodel == model
        return crit_pure(model,nothing)
    else
        return crit_pure(satmodel) #overloaded somewhere
    end
end

function crit_x_to_v(lbv,x)
    lo = 0.001
    hi = 1.0
    η = lo + (hi - lo)/(exp(-x) + 1)
    return lbv/η
end

function crit_v_to_x(lbv,v)
    η = lbv/v
    lo = 0.001
    hi = 1.0
    return -log((hi - lo)/(η - lo) - 1)
end

function crit_pure(model::EoSModel,x0,z = SA[1.0];options = NEqOptions())
    return critical_point_solver(model,x0,z;options)
end

function critical_point_solver(model::EoSModel,x0,z = SA[1.0];options = NEqOptions())
    check_arraysize(model,z)
    #f! = (F,x) -> obj_crit(model, F, x[1]*T̄, exp10(x[2]))
    zp = primalval(z)
    primalmodel = primalval(model)
    if x0 === nothing
        _x0 = x0_crit_pure(primalmodel,zp)
    elseif x0 isa Type && x0 <: Number
        x0x = x0_crit_pure(primalmodel,zp)
        _x0 = x0.(x0x)
    else
        _x0 = x0
    end
    x01,x02 = _x0
    T̄  = T_scale(primalmodel,zp)*one(x01*one(x02))
    Tc0 = x01*T̄
    lbv0 = lb_volume(primalmodel,Tc0,zp)
    _1 = oneunit(T̄)
    Tx0 = primalval(x01)
    vc0 = exp10(primalval(x02))
    x0 = vec2(Tx0,crit_v_to_x(lbv0,vc0),_1)
    zz = zp/sum(zp)
    f! = WithContext(CritPure(primalmodel,primalval(T̄),zz),∂Tag{critical_point_solver}())
    solver_res = Solvers.nlsolve(f!, x0, TrustRegion(Newton(), NLSolvers.NWI()), options)
    r  = Solvers.x_sol(solver_res)
    !__check_convergence(solver_res) && (r .= NaN)
    T_c = r[1]*T̄
    lbv = lb_volume(primalmodel,T_c,zp)
    V_c = crit_x_to_v(lbv,r[2])
    p_c = pressure(primalmodel, V_c, T_c, zp)
    if p_c < 0
        p_c *= NaN
        V_c *= NaN
        T_c *= NaN
    end
    crit = (T_c, p_c, V_c)
    tup = (model,z)
    λtup = (primalmodel,zp)
    return crit_pure_ad(crit,tup,λtup)
end

function crit_pure_ad(crit,tup,λtup)
    if any(has_dual,tup) # do check here to avoid recomputation of pressure if no AD
        λx = SVector(crit[1],crit[3])
        f(x,tups) = begin
            model,z = tups
            Tc,Vc = x
            _,F1,F2 = p∂p∂2p(model, Vc, Tc, z)
            return SVector(F1,F2)
        end
        T_c,V_c = __gradients_for_root_finders(λx,tup,λtup,f)
        P_c = pressure(model,V_c,T_c,z)
        return (T_c,P_c,V_c)
    else
        return crit
    end
end

struct CritPure{M,T,Z}
    model::M
    T̄::T
    z::Z
end

StaticForwardDiffTags.deferred_valtype(m::CritPure) = Base.promote_eltype(m.model,m.T̄,m.z)

function ObjCritPure(model::T,F,T̄,x,z) where T
    sv = __ObjCritPure(model,T̄,x,z)
    F[1] = sv[1]
    F[2] = sv[2]
    return F
end

function (crit::CritPure{M,T,Z})(F,x) where {M,T,Z}
    sv = __ObjCritPure(crit.model,crit.T̄,x,crit.z)
    F[1] = sv[1]
    F[2] = sv[2]
    return F
end

function (crit::CritPure{M,T,Z})(x) where {M,T,Z}
    return __ObjCritPure(crit.model,crit.T̄,x,crit.z)
end

function __ObjCritPure(model::T,T̄,x,z) where T
    T_c = x[1]*T̄
    lbv = lb_volume(model,T_c,z)
    V_c = crit_x_to_v(lbv,x[2])
    RT = Rgas(model)*T_c
    ∂²A∂V²_scale = V_c*V_c/RT
    ∂³A∂V³_scale = ∂²A∂V²_scale*V_c
    _, ∂²A∂V², ∂³A∂V³ = p∂p∂2p(model, V_c, T_c, z)
    F1 = -∂²A∂V²*∂²A∂V²_scale
    F2 = -∂³A∂V³*∂³A∂V³_scale
    return SVector((F1,F2))
end

"""
    mechanical_critical_point(model,z = SA[1.0],x0 = x0_crit_pure(model,z))

Returns the mechanical critical point of an `EoSModel`. The mechanical critical point is the point where the first and second derivatives of the pressure function vanish.
For single component models, the function is equivalent to [`crit_pure`](@ref).
Returns a tuple, containing:

- Mechanical critical temperature `[K]`
- Mechanical critical pressure `[Pa]`
- Mechanical critical molar volume `[m³·mol⁻¹]`
"""
function mechanical_critical_point end

mechanical_critical_point(model,z,x0) = crit_pure(model,x0,z)
mechanical_critical_point(model,z) = mechanical_critical_point(model,z,x0_crit_pure(model,z))
mechanical_critical_point(model) = crit_pure(model)

"""
    mechanical_crit_bracket(model, z=SA[1.0], T0 = T_scale(model,z); tol = 1e-2, max_iters = 60, max_widen = 5)



Find a temperature bracket `(Tmin, Tmax)` such that a coexistence loop exists at `Tmin` but not at `Tmax`, starting from the probe temperature `T0`.

Convergence is declared when the temperature bracket has collapsed (`tol`, relative on `τ = 1/T`).

## keyword Inputs

- `tol`             : Relative tolerance on `τ = 1/T`.
- `max_iters`       : Maximum iterations of the algorithm. If convergence is not reached but there is a valid bracket, then the function will return the tightest bracket calculated.
- `max_widen`       : if `T0` does not contain a coexistence loop, then `Clapeyron.probe_subcritical` is called with `max_iters == max_widen`.

## Outputs

 - `(Tmin, Tmax)`   : Final temperature bracket. Coexistence exists at `Tmin`
                      (so `Tmin < Tc`) and not at `Tmax` (so `Tmax > Tc`).
 - `(vmin, vmax)`   : Bracket for the critical volume, taken from the stable
                      coexistence volumes at `Tmin` (the last
                      temperature where a loop was found). Because `vmin` and
                      `vmax` both converge to `vc` as `T → Tc`, this is a
                      legitimate bracket of the critical volume.
 - `vbest`          : Best volume estimate. unstable middle volume at `T = Tmin`.
"""
function mechanical_crit_bracket(model::EoSModel, z=SA[1.0], T0 = T_scale(model, z);
                                tol_T::Float64      = 1e-2,
                                max_iters::Int      = 60,
                                max_widen::Int      = 60)



    # --- step 0: initial temperature bracket --------------------------
    # widen Tlo downwards until a Maxwell loop appears

    sl0,sv0,sm0,status0 = coexistence_volumes_bracketed(model,T0,z)
    TT = Base.promote_eltype(model, T0, z ,sl0.v, sv0.v,sm0.v)

    if status0 === :failure
        T0_bounds, sl0, sv0, sm0 = probe_subcritical(model, sm0.v, T0, z, max_iters = max_widen) #solver to get a subcritical point
        Tlo,Thi = TT(T0_bounds[1]),TT(T0_bounds[2])
    else
        Tlo = TT(T0)
        Thi = TT(Inf)
    end

    sm_old = (v = TT(NaN),p = TT(NaN),dpdv = TT(NaN),d2pdv2 = TT(NaN))
    Told = TT(NaN)
    sl,sv,sm = sl0,sv0,sm0
    Tm = TT(NaN)

    for i in 1:max_iters

        #generate new temperature estimate (guarded secant on dpdv)
        τlo,τhi = 1/Tlo,1/Thi
        τm = 0.5*(τlo + τhi)
        Tmid = 1 / τm
        if !isnan(sm_old.v) && isfinite(Thi)
            τ_secant = τlo - (τlo - 1/Told)/(sm.dpdv - sm_old.dpdv)
            T_secant = 1/τ_secant
            if Tlo <= T_secant <= Thi
                #damp secant so it does not get stuck in a subcritical spinodal
                Tm = 0.9/τ_secant + 0.1*Tmid
            else
                Tm = Tmid
            end
        else
            Tm = Tmid
        end

        sli, svi, smi, status = coexistence_volumes_bracketed(model, Tm, z, sl.v, sv.v, sm.v)

        lo_updated = false

        if status !== :failure && smi.dpdv > 0
            sm_old,Told = sm,Tlo
            sl,sv,sm = sli,svi,smi
            Tlo = Tm
            lo_updated = true
        else
            Thi = Tm
        end

        if lo_updated
            #contract the high temperature bound
            Trange = range(Thi,Tlo,3)
            for k in 1:4
                Tm = Trange[k+1]
                slk, svk, smk, status = coexistence_volumes_bracketed(model, Tm, z, sl.v, sv.v, sm.v)
                if status !== :failure
                    sm_old,Told = sm,Tlo
                    sl,sv,sm = slk,svk,smk
                    Tlo = Tm
                    break
                else
                    #high bound contracted
                    Thi = Tm
                end
            end
        end
        #@show Tlo,Thi
        conv_T, _ = Solvers.convergence(1/Tlo, 1/Thi, zero(τm), tol_T)
        conv_T && break
    end

    if !isfinite(Tlo) || !isfinite(Thi)
        throw(ArgumentError(lazy"Clapeyron.mechanical_crit_bracket: no successful coexistence bracket was ever produced."))
    end

    return (Tlo, Thi), (sl.v, sv.v), sm.v
end

"""
    x0_crit_bisection(model::EoSModel, z=SA[1.0], T0 = nothing; kwargs...)

Returns a initial point for [`crit_pure`](@ref) and [`crit_mix`](@ref), by using a bracketing algorithm ([`mechanical_crit_bracket`](@ref)) to locate the mechanical critical point.
Useful for models with hard-to converge critical points, or for models where even an approximation of the mechanical critical coordinates are not known beforehand.

Used by default for non-splittable models and electrolyte models as initial point in `crit_mix`.

`T0` can be used to provide an initial point for the bisection procedure.

Using this method as a default for a defined `EoSModel` just requires adding a method to `x0_crit_pure` (or `x0_crit_mix`)

```julia
Clapeyron.x0_crit_pure(model::MyModel,z) = Clapeyron.x0_crit_bisection(model,z)
Clapeyron.x0_crit_mix(model::MyModel,z)  = Clapeyron.x0_crit_bisection(model,z)
```
"""
function x0_crit_bisection(model::EoSModel, z=SA[1.0], T0 = nothing; kwargs...)
    Ts = T_scale(model, z)
    if T0 === nothing
        TT, _, v_best = mechanical_crit_bracket(model, z; kwargs...)
    else
        TT, _, v_best = mechanical_crit_bracket(model, z, T0; kwargs...)
    end
    Tlo, Thi = TT
    T0 = 0.5*(Tlo + Thi)
    return (T0 / Ts, log10(v_best))
end

function probe_subcritical(model::EoSModel, v0, T0, z = SA[1.0]; max_iters::Int = 40)
    TT = Base.promote_eltype(model, v0, T0, z)
    T  = TT(T0)
    v  = TT(v0)
    R  = Rgas(model)
    g  = TT(NaN)
    δlo, δhi = TT(1.01748), TT(2.0)        # static bounds on the tau growth factor
    Tlo, Thi = T, TT(Inf)              # Inf: a loop at T0 means T0 is a valid lower bound

    for i in 1:max_iters
        lb = lb_volume(model, Tlo, z)
        # confirmation + seeds for the main loop
        sl, sv, sm, status = coexistence_volumes_bracketed(model, Tlo, z, sqrt(lb*v), 10*v, v)
        if status !== :failure && sm.dpdv > 0
            return (Tlo, Thi), sl, sv, sm
        end
        Thi = Tlo

        # hill-climb v toward the best probe (ignore NaN dpdv)
        sx = (sl, sv, sm)
        g_new, idx = findmax(w -> isnan(w.dpdv) ? TT(-Inf) : w.dpdv, sx)
        if isfinite(g_new) && (isnan(g) || g_new > g)
            v = sx[idx].v
            g = g_new
        end

        # descend T: translated-ideal-gas proxy, with the step bounded in tau
        τ = 1/Tlo
        if isnan(g)
            τnew = δlo*τ                       # no usable g yet: forced descent
        else
            ΔT       = g*(v - lb)^2/R          # g <= 0 here, so ΔT <= 0
            τ_update = 1/(Tlo + ΔT)            # negative if Tlo + ΔT <= 0: falls in the floor branch
            τnew     = clamp(τ_update, δlo*τ, δhi*τ)
        end
        Tlo = 1/τnew
    end

    throw(ArgumentError(lazy"probe_subcritical: no unstable point found in $max_iters iterations (last T = $Tlo)"))
end


@public x0_crit_bisection,mechanical_crit_bracket