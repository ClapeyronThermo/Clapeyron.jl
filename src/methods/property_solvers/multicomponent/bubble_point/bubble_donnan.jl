#=
Donnan-equilibrium phase-equilibrium solver, for any `ElectrolyteModel`
whose components carry a well-defined net charge (via `component_charges`).
Given `T` and one phase's full, electroneutral composition `x`, solves
jointly for the other (coexisting) phase's full composition `y`, both
phases' volumes, and the Donnan potential `Ψ`, via `N+2` equations:
  - 1 mechanical equilibrium:  P(V_dense,x) = P(V_dilute,y)
  - N electrochemical potential matching (one per species):
        μᵢ(V_dense,x) = μᵢ(V_dilute,y) + Zᵢ·Rgas·T·Ψ
  - 1 electroneutrality of the *solved* phase: Σ Zᵢ yᵢ = 0
against `N+2` unknowns: V_dense, V_dilute, the dilute phase's `N-1`
independent mole fractions, and Ψ.

`bubble_pressure`'s generic `EoSModel`-typed dispatch chain
(`index_reduction`, `bubble_pressure_ad`) assumes a standard `(P,vl,vv,y)`
return and has no room for the extra `Ψ` unknown or the extra
electroneutrality equation, so `DonnanBubblePressure`'s own `bubble_pressure`
method dispatches directly on `model::ElectrolyteModel` -- more specific
than, and so preferred over, the generic `model::EoSModel` method -- rather
than going through `bubble_pressure_impl`.
=#

"""
    DonnanBubblePressure(; u0=nothing, xtol=1e-10, max_iters=1000, verbose=false, psi_zero=false)

Method type for `bubble_pressure(model::ElectrolyteModel, T, x, ::DonnanBubblePressure)`.
See this file's header for the full derivation.

`u0`, if given, is the full initial guess -- `[V_dense, V_dilute, y_free...,
Ψ]` (`N+2` entries) normally, or `[V_dense, V_dilute, y_free...]` (`N+1`
entries, no `Ψ`) when `psi_zero=true`; if `nothing`, a guess is constructed
from a packing-fraction bracket.

`psi_zero=true` selects the reduced equation set for systems where `Ψ=0`
is known a priori -- e.g. a net-neutral polymer/solvent species plus a
symmetric 1:1 salt (no counter-ion needed for the neutral species' own
electroneutrality), where the model's cation↔anion charge-swap symmetry
forces `Ψ≡0` exactly. With `Ψ` fixed at 0 rather than solved for, the
electrochemical-potential-matching equations reduce to ordinary chemical-
potential matching, and the problem becomes mathematically identical to a
standard (uncharged) N-component bubble-pressure calculation.
"""
struct DonnanBubblePressure <: BubblePointMethod
    u0::Union{Nothing,Vector{Float64}}
    xtol::Float64
    max_iters::Int
    verbose::Bool
    psi_zero::Bool
end

function DonnanBubblePressure(; u0=nothing, xtol::Float64=1e-10, max_iters::Int=1000, verbose::Bool=false, psi_zero::Bool=false)
    return DonnanBubblePressure(u0, xtol, max_iters, verbose, psi_zero)
end

_donnan_mole_fractions_from_free(free::AbstractVector) = vcat(free, 1 - sum(free))

function donnan_bubble_residual!(F::AbstractVector, u::AbstractVector, model::ElectrolyteModel, T, x_dense::AbstractVector, Z::AbstractVector{<:Integer})
    N = length(Z)
    Vd, Vl = u[1], u[2]
    y = _donnan_mole_fractions_from_free(view(u, 3:1+N))
    Ψ = u[end]

    F[1] = pressure(model, Vd, T, x_dense) - pressure(model, Vl, T, y)
    μd = VT_chemical_potential(model, Vd, T, x_dense)
    μl = VT_chemical_potential(model, Vl, T, y)
    RT = Rgas(model) * T
    for i in 1:N
        F[1+i] = (μd[i] - μl[i]) / RT - Z[i] * Ψ
    end
    F[N+2] = dot(Z, y)
    return F
end

function donnan_bubble_residual_psi_zero!(F::AbstractVector, u::AbstractVector, model::ElectrolyteModel, T, x_dense::AbstractVector, N::Int)
    Vd, Vl = u[1], u[2]
    y = _donnan_mole_fractions_from_free(view(u, 3:1+N))

    F[1] = pressure(model, Vd, T, x_dense) - pressure(model, Vl, T, y)
    μd = VT_chemical_potential(model, Vd, T, x_dense)
    μl = VT_chemical_potential(model, Vl, T, y)
    RT = Rgas(model) * T
    for i in 1:N
        F[1+i] = (μd[i] - μl[i]) / RT
    end
    return F
end

function _donnan_bracket_coexistence_volumes(model::ElectrolyteModel, T::Real, x::AbstractVector; nV::Int=400, Vmax_factor::Float64=1e6)
    lbv = lb_volume(model, x)
    Vs = exp.(range(log(lbv * 1.01), log(lbv * Vmax_factor), length=nV))
    Ps = [pressure(model, V, T, x) for V in Vs]
    i = 2
    while i <= nV && Ps[i] < Ps[i-1]
        i += 1
    end
    i > nV && throw(ErrorException("pressure strictly decreasing over the whole scanned range: no loop found at T=$T for this composition (not below the fixed-ray critical point, or Vmax_factor too small)"))
    i_min = i - 1
    while i <= nV && Ps[i] > Ps[i-1]
        i += 1
    end
    i > nV && throw(ErrorException("pressure never resumed decreasing after the loop's rise at T=$T: scan range too narrow, increase Vmax_factor"))
    i_max = i - 1
    return Vs[i_min] * 0.5, Vs[i_max] * 2.0
end

function _default_donnan_u0(model::ElectrolyteModel, T::Real, x_dense::AbstractVector, Z::AbstractVector{<:Integer})
    N = length(Z)
    Vd0, Vl0 = _donnan_bracket_coexistence_volumes(model, T, x_dense)
    return vcat(Vd0, Vl0, x_dense[1:N-1], 0.0)
end

"""
    donnan_psi_bulk(model::ElectrolyteModel, T, ρbulk_dense::AbstractVector, ρbulk_dilute::AbstractVector)

Read off the Donnan potential `Ψ` (dense-minus-dilute convention, matching
`donnan_bubble_residual!`'s own `μᵢ(dense) - μᵢ(dilute) = Zᵢ·Rgas·T·Ψ`)
between two ALREADY-KNOWN coexisting bulk states, given directly as
per-component densities rather than `(V,x)` pairs -- since chemical
potential only depends on density (`VT_chemical_potential(model,1,T,ρ)` and
`VT_chemical_potential(model,V,T,ρ*V)` agree for any `V`, confirmed
directly), no volume/mole-fraction bookkeeping is needed from the caller.

No iterative solve: unlike `bubble_pressure(...,DonnanBubblePressure())`,
which solves the FULL coupled equilibrium (`Ψ` along with both phases'
volumes/composition), this is a direct readout for a caller that already
has a converged coexisting pair of bulk states in hand (e.g. a DFT solver's
own `structure.ρbulk`/`structure.topology.ρbulk2`, built earlier from a
`bubble_pressure` result) and only needs `Ψ` itself, not to re-solve the
equilibrium that produced it. Averaged over every charged component
(`Zᵢ≠0`) for numerical robustness -- at a genuine coexistence point these
should all agree; returns exactly `0.0` if `model` has no charged
components (nothing to average).
"""
function donnan_psi_bulk(model::ElectrolyteModel, T, ρbulk_dense::AbstractVector, ρbulk_dilute::AbstractVector)
    Z = component_charges(model)
    charged = findall(!iszero, Z)
    isempty(charged) && return 0.0
    μd = VT_chemical_potential(model, 1.0, T, ρbulk_dense)
    μl = VT_chemical_potential(model, 1.0, T, ρbulk_dilute)
    RT = Rgas(model) * T
    return sum(i -> (μd[i] - μl[i]) / (RT * Z[i]), charged) / length(charged)
end

"""
    bubble_pressure(model::ElectrolyteModel, T, x, method::DonnanBubblePressure=DonnanBubblePressure())

Donnan-equilibrium-augmented bubble pressure: given `T` and one phase's
full electroneutral composition `x`, returns `(P, V_dense, V_dilute, y, Ψ)`.
With `method.psi_zero=true`, solves the reduced `Ψ=0` system instead (see
`DonnanBubblePressure`'s docstring).
"""
function bubble_pressure(model::ElectrolyteModel, T, x::AbstractVector, method::DonnanBubblePressure=DonnanBubblePressure())
    Z = component_charges(model)
    N = length(Z)
    length(x) == N || throw(ArgumentError("x must have length $N (one entry per component), got $(length(x))"))
    isapprox(sum(x), 1.0; atol=1e-8) || throw(ArgumentError("x must sum to 1 (mole fractions), got sum=$(sum(x))"))
    isapprox(dot(Z, x), 0.0; atol=1e-6) || throw(ArgumentError("x must be electroneutral (Σ Zᵢxᵢ=0); got $(dot(Z, x))"))

    if method.psi_zero
        u0_full = method.u0 === nothing ? _default_donnan_u0(model, T, x, Z) : copy(method.u0)
        u0 = length(u0_full) == N + 1 ? u0_full : u0_full[1:N+1]  # drop Ψ if a full-length guess was supplied
        f_psi0!(F, u) = donnan_bubble_residual_psi_zero!(F, u, model, T, x, N)
        result = Solvers.nlsolve(f_psi0!, u0)
        u = Solvers.x_sol(result)

        converged = __check_convergence(result) && !any(isnan, u)
        converged || (u .= NaN)

        Vd, Vl = u[1], u[2]
        y = _donnan_mole_fractions_from_free(view(u, 3:1+N))
        if converged && !isapprox(dot(Z, y), 0.0; atol=1e-6)
            @warn "DonnanBubblePressure(psi_zero=true): converged y is not electroneutral (Σ Zᵢyᵢ=$(dot(Z, y))) -- Ψ=0 was likely not a valid assumption for this model/composition; use psi_zero=false."
        end
        P = converged ? pressure(model, Vd, T, x) : NaN
        return (P=P, V_dense=Vd, V_dilute=Vl, y=y, Ψ=0.0)
    end

    u0 = method.u0 === nothing ? _default_donnan_u0(model, T, x, Z) : copy(method.u0)
    f!(F, u) = donnan_bubble_residual!(F, u, model, T, x, Z)
    result = Solvers.nlsolve(f!, u0)
    u = Solvers.x_sol(result)

    # At a species' mole fraction of exactly 0, the ideal mixing entropy's
    # log(xᵢ) term diverges, propagating NaN through the whole residual --
    # nlsolve can silently "converge" to a NaN-contaminated point, so check
    # explicitly rather than trust whatever x_sol returns.
    converged = __check_convergence(result) && !any(isnan, u)
    converged || (u .= NaN)

    Vd, Vl, Ψ = u[1], u[2], u[end]
    y = _donnan_mole_fractions_from_free(view(u, 3:1+N))
    P = converged ? pressure(model, Vd, T, x) : NaN
    return (P=P, V_dense=Vd, V_dilute=Vl, y=y, Ψ=Ψ)
end

"""
    find_critical_salt_fraction(model::EoSModel, T; s_lo=0.01, s_hi=0.999, max_iters=60)

Critical salt fraction `s_c` for a net-neutral species + symmetric 1:1 salt
system (`x(s) = [1-s, s/2, s/2]`, `s` the total salt mole fraction) at a
FIXED temperature `T` -- the composition-space analogue of a critical
*temperature* at fixed composition, with the roles of `T` and composition
swapped: here `T` is held fixed and the composition ray `x(s)` is what's
bisected on.

Bisects `s` on whether `crit_mix(model,x(s))`'s own critical temperature
`Tc_mix(s)` is above or below the target `T` -- `crit_mix` solves the full
multicomponent critical-point condition (stability against both density and
composition perturbations), which is the physically-relevant criterion for
whether `x(s)` itself sits at a genuine multicomponent critical composition
at `T`. `Tc_mix(s) > T` means still unstable (coexistence exists somewhere
reachable from `x(s)` at this `T`); `Tc_mix(s) < T` means stable. `crit_mix`
becomes numerically fragile very close to `s→1` (pure salt, a genuinely
degenerate limit) -- a `NaN` return during bisection is treated as "still
unstable" (the fragility is on the far side of where `s_c` is expected for
any solvent/polymer-containing system).

Returns `(s_c, V_c, P_c)`.
"""
function find_critical_salt_fraction(model::EoSModel, T::Real; s_lo::Float64=0.01, s_hi::Float64=0.999, max_iters::Int=60)
    xs(s) = [1 - s, s / 2, s / 2]

    crit_mix_at(s) = crit_mix(model, xs(s))

    Tc_lo, = crit_mix_at(s_lo)
    Tc_hi, = crit_mix_at(s_hi)
    isnan(Tc_lo) || Tc_lo > T || throw(ArgumentError("no coexistence found even at s_lo=$s_lo (crit_mix Tc=$Tc_lo < T=$T); lower it"))
    isnan(Tc_hi) || Tc_hi < T || throw(ArgumentError("still coexisting at s_hi=$s_hi (crit_mix Tc=$Tc_hi > T=$T); raise it"))

    lo, hi = s_lo, s_hi
    sc, Pc, Vc = NaN, NaN, NaN
    for _ in 1:max_iters
        mid = 0.5 * (lo + hi)
        Tc_mid, Pc_mid, Vc_mid = crit_mix_at(mid)
        if isnan(Tc_mid) || Tc_mid < T
            hi = mid
        else
            lo = mid
        end
        sc, Pc, Vc = mid, Pc_mid, Vc_mid
    end
    return (s_c=sc, V_c=Vc, P_c=Pc)
end

export DonnanBubblePressure, donnan_bubble_residual!, donnan_bubble_residual_psi_zero!
export find_critical_salt_fraction, donnan_psi_bulk
