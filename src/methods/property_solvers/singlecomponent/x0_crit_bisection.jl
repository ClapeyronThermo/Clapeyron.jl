"""
    has_physical_loop(model::EoSModel, T, z=SA[1.0]; nV=400, Vmax_factor=1e5)

Scan pressure over the physical volume range (`V > lb_volume`, i.e. packing
fraction `η<1`) at fixed `T` along composition ray `z`, and report whether
`P(V)` is non-monotonic there -- the direct, derivative-free signature of a
genuine two-phase region. Ignores any apparent non-monotonicity below
`lb_volume`, which is just a repulsive term diverging at unphysical
densities, not a real phase transition.
"""
function has_physical_loop(model::EoSModel, T::Real, z=SA[1.0]; nV::Int=400, Vmax_factor::Float64=1e5)
    lbv = lb_volume(model, T, z)
    Vs = exp.(range(log(lbv * 1.01), log(lbv * Vmax_factor), length=nV))
    prev = pressure(model, Vs[1], T, z)
    for V in Vs[2:end]
        p = pressure(model, V, T, z)
        p > prev && return true
        prev = p
    end
    return false
end

"""
    x0_crit_bisection_TV(model::EoSModel, z=SA[1.0]; Tlo_factor=0.1, Thi_factor=10.0, max_widen=8, eta0=0.02)

Bisect on `T` (geometric) for the onset of a genuine physical pressure loop
along the fixed composition ray `z` (via [`has_physical_loop`](@ref)), and
pair it with a low-packing-fraction volume guess `Vc0 = lb_volume(model,Tc0,z)/eta0`.

Intended for models whose critical packing fraction is far from the
~0.13-0.3 range the generic `x0_crit_pure`/`x0_crit_mix` defaults assume
(e.g. electrostatically-driven phase separation, critical packing fraction
~0.005-0.02) -- call this from a model-specific `x0_crit_pure`/`x0_crit_mix`
override (see [`x0_crit_pure_bisection`](@ref)/[`x0_crit_mix_bisection`](@ref))
rather than reaching for a bespoke solver wrapper outside Clapeyron.

Returns `(Tc0, Vc0)`: a warm-start pair for `crit_pure`/`crit_mix`, not the
converged critical point itself.
"""
function x0_crit_bisection_TV(model::EoSModel, z=SA[1.0]; Tlo_factor::Float64=0.1, Thi_factor::Float64=10.0,
                               max_widen::Int=8, eta0::Float64=0.02)
    Ts = T_scale(model, z)
    Thi = Ts * Thi_factor
    Tlo = Ts * Tlo_factor
    has_physical_loop(model, Thi, z) && throw(ArgumentError("no single-phase region found at Thi=$Thi; widen Thi_factor"))
    for _ in 1:max_widen
        has_physical_loop(model, Tlo, z) && break
        Tlo /= 3
    end
    has_physical_loop(model, Tlo, z) || throw(ArgumentError("no coexistence found even at Tlo=$Tlo; widen Tlo_factor or max_widen"))

    for _ in 1:60
        Tmid = sqrt(Thi * Tlo) # geometric bisection (T spans orders of magnitude)
        if has_physical_loop(model, Tmid, z)
            Tlo = Tmid
        else
            Thi = Tmid
        end
    end
    Tc0 = sqrt(Thi * Tlo)
    Vc0 = lb_volume(model, Tc0, z) / eta0
    return Tc0, Vc0
end

"""
    x0_crit_pure_bisection(model::EoSModel, z=SA[1.0]; kwargs...)

`x0_crit_pure`-shaped warm start (`(Tc0/T_scale(model,z), log10(Vc0))`),
built from [`x0_crit_bisection_TV`](@ref). A model whose critical point the
generic `x0_crit_pure` default can't find (a spurious stationary point from
its own ~0.3-packing-fraction guess) should define
`x0_crit_pure(model::MyModel,z) = x0_crit_pure_bisection(model,z)`.
"""
function x0_crit_pure_bisection(model::EoSModel, z=SA[1.0]; kwargs...)
    Ts = T_scale(model, z)
    Tc0, Vc0 = x0_crit_bisection_TV(model, z; kwargs...)
    return (Tc0 / Ts, log10(Vc0))
end

export has_physical_loop, x0_crit_bisection_TV, x0_crit_pure_bisection
