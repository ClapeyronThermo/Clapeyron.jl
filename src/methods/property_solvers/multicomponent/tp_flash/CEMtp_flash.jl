"""
    CEMFlash(;δ = nothing, equilibrium = :auto, polish = MichelsenTPFlash(), threaded = true, diagram = nothing)

Method to solve the non-reactive multicomponent, multiphase flash problem with the convex envelope method (CEM) [1].

The Gibbs energy of every aggregate state is evaluated on a lattice of the composition simplex (step `1/δ`). The phase split is the facet of the lower convex envelope above the feed composition, found with a linear program (see [`CEMPhaseDiagram`](@ref)).
The number of phases does not need to be known, and the method does not depend on initial guesses.
The CEM result has a discretization error of order `1/δ`. By default, it is refined with a local flash method (`polish`).

If a [`Clapeyron.CEMDiagram`](@ref) obtained with [`phase_diagram`](@ref) is passed via `diagram`, the stored envelope is used and no Gibbs energy evaluation on the grid is necessary (the conditions and components must match).

## Keyword arguments
- `δ`: number of divisions of the composition axes. By default, `256`, `128` and `64` for binary, ternary and quaternary mixtures (as in [1]), and `32` for more components. Ignored if `diagram` is given.
- `equilibrium`: `:auto` considers liquid and vapour phases, `:lle` only liquid phases. Ignored if `diagram` is given.
- `polish`: `TPFlashMethod` used to refine the CEM phases (initialized with the CEM result), or `nothing` to return the CEM result. The refined result is only accepted if it has a lower Gibbs energy than the CEM result. If `MichelsenTPFlash` is given and the CEM finds more than two phases, a `MultiPhaseTPFlash` with the same tolerances is used. Results with solid phases are not refined.
- `threaded`: evaluate the grid in parallel.
- `diagram`: a precomputed `CEMDiagram`.

Pure solid phases are considered for `CompositeModel`s with a `SolidHfus` solid model. Solid phases are returned with a `NaN` volume.

## References
1. Göttl, Q., Rosen, N., & Burger, J. (2025). Convex envelope method for T, p flash calculations for mixtures with an arbitrary number of components and arbitrary aggregate states. Computers & Chemical Engineering, 203, 109326. [doi:10.1016/j.compchemeng.2025.109326](https://doi.org/10.1016/j.compchemeng.2025.109326)
"""
Base.@kwdef struct CEMFlash <: TPFlashMethod
    δ::Union{Int,Nothing} = nothing
    equilibrium::Symbol = :auto
    polish::Union{Nothing,TPFlashMethod} = MichelsenTPFlash()
    threaded::Bool = true
    diagram::Union{Nothing,CEMDiagram} = nothing
end

#the stored envelope contains all components
supports_reduction(method::CEMFlash) = method.diagram === nothing

function index_reduction(method::CEMFlash,idx::AbstractVector)
    polish = method.polish === nothing ? nothing : index_reduction(method.polish,idx)
    return CEMFlash(method.δ,method.equilibrium,polish,method.threaded,method.diagram)
end

function tp_flash_impl(model::EoSModel,p,T,z,method::CEMFlash)
    N = length(model)
    d = method.diagram
    if d === nothing
        equilibrium = method.equilibrium
        env = cem_envelope(model,p,T,something(method.δ,cem_default_δ(N)),equilibrium,method.threaded)
        stable = CEMStability(env)
    else
        component_list(model) == d.components || throw(ArgumentError("CEMFlash: the components of the model and of the diagram do not match"))
        (isapprox(p,d.p) && isapprox(T,d.T)) || throw(ArgumentError("CEMFlash: the diagram was calculated at p = $(d.p), T = $(d.T), but the flash is at p = $p, T = $T"))
        equilibrium = d.equilibrium
        env = d.envelope
        stable = CEMStability(env,d.point_facet)
    end
    fluidmodel,fluid,solid = cem_models(model,p,T,z,equilibrium)
    ws = CEMWorkspace(env)
    zf = collect(Float64,z)
    cem_reset!(ws,env,zf)
    cem_lp!(ws,env,zf)
    result = cem_flash_result(env,fluid,p,T,zf,ws,stable)
    return cem_polish(fluidmodel,p,T,zf,result,method.polish)
end

cem_polish(model,p,T,z,result,::Nothing) = result

function cem_polish(model,p,T,z,result,method::TPFlashMethod)
    #single phase, or solid phases present (the local methods do not handle solids)
    (numphases(result) < 2 || !all(isfinite,result.volumes)) && return result
    m = cem_polish_method(method,result)
    m === nothing && return result
    polished = try
        tp_flash_impl(model,p,T,z,m)
    catch
        return result
    end
    all(isfinite,polished.volumes) && all(x -> all(isfinite,x),polished.compositions) || return result
    #reduced molar gibbs energies (same reference). The discrete envelope lies above the exact one: a correct refinement can only lower g
    g0,g1 = result.data.g,polished.data.g
    return isfinite(g1) && g1 <= g0 + 1e-10*(1 + abs(g0)) ? polished : result
end

cem_polish_method(method::TPFlashMethod,result) = method

function cem_polish_method(m::MichelsenTPFlash,result)
    if numphases(result) > 2
        return MultiPhaseTPFlash(;K_tol = m.K_tol,ss_iters = m.ss_iters,nacc = m.nacc,second_order = m.second_order,verbose = m.verbose,flash_result = result)
    end
    equilibrium = m.equilibrium
    if is_unknown(equilibrium)
        equilibrium = result.data.vapour_idx > 0 ? :vle : :lle
    end
    return MichelsenTPFlash(;equilibrium,K_tol = m.K_tol,ss_iters = m.ss_iters,nacc = m.nacc,second_order = m.second_order,
        noncondensables = m.noncondensables,nonvolatiles = m.nonvolatiles,verbose = m.verbose,flash_result = result)
end

function cem_polish_method(m::MultiPhaseTPFlash,result)
    return MultiPhaseTPFlash(;K_tol = m.K_tol,ss_iters = m.ss_iters,nacc = m.nacc,second_order = m.second_order,
        max_phases = m.max_phases,phase_iters = m.phase_iters,verbose = m.verbose,flash_result = result)
end

export CEMFlash
