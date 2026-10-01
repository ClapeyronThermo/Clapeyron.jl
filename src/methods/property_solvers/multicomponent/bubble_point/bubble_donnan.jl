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
    y = FractionVector(view(u, 3:1+N))
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
    y = FractionVector(view(u, 3:1+N))

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

function bubble_pressure_impl(model::ElectrolyteModel, T, nn::AbstractVector, method::DonnanBubblePressure)
    x = nn ./ sum(nn)
    Z = component_charges(model)
    N = length(Z)
    length(x) == N || throw(ArgumentError("x must have length $N (one entry per component), got $(length(x))"))

    electroneutral_check(@view(ν[i,:]),charges)
    nan = Base.promote_eltype(model,T,x)
    if method.psi_zero
        u0_full = method.u0 === nothing ? _default_donnan_u0(model, T, x, Z) : copy(method.u0)
        u0 = length(u0_full) == N + 1 ? u0_full : u0_full[1:N+1]  # drop Ψ if a full-length guess was supplied
        f_psi0!(F, u) = donnan_bubble_residual_psi_zero!(F, u, model, T, x, N)
        result = Solvers.nlsolve(f_psi0!, u0)
        u = Solvers.x_sol(result)

        converged = __check_convergence(result) && !any(isnan, u)
        converged || (u .= NaN)

        Vd, Vl = u[1], u[2]
        yn = view(u, 3:1+N)
        y = vcat(yn,1 - sum(yn))
        if converged && !isapprox(dot(Z, y), 0.0; atol=1e-6)
            @warn "DonnanBubblePressure(psi_zero=true): converged y is not electroneutral (Σ Zᵢyᵢ=$(dot(Z, y))) -- Ψ=0 was likely not a valid assumption for this model/composition; use psi_zero=false."
        end
        P = converged ? pressure(model, Vd, T, x) : nan
        return P,Vl,Vd,y
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
    yn = view(u, 3:1+N)
    y = vcat(yn,1 - sum(yn))
    P = converged ? pressure(model, Vd, T, x) : nan
    return P,Vd,Vl,y
end

export DonnanBubblePressure, donnan_bubble_residual!, donnan_bubble_residual_psi_zero!
