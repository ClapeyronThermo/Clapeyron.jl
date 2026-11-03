"""
    CEMFlash(;δ = nothing, equilibrium = :auto, polish = MichelsenTPFlash(second_order = true), threaded = true, diagram = nothing)

Method to solve the non-reactive multicomponent, multiphase flash problem with the convex envelope method (CEM) [1].

The Gibbs energy of every aggregate state is evaluated on a lattice of the composition simplex (step `1/δ`). The phase split is the facet of the lower convex envelope above the feed composition, found with a linear program (see [`CEMPhaseDiagram`](@ref)).
The number of phases does not need to be known, and the method does not depend on initial guesses.
The CEM result has a discretization error of order `1/δ`. By default, it is refined with a local flash method (`polish`).

If a `CEMDiagram` obtained with [`phase_diagram`](@ref) is passed via `diagram`, the stored envelope is used and no Gibbs energy evaluation on the grid is necessary (the conditions and components must match).

## Keyword arguments
- `δ`: number of divisions of the composition axes. By default, `256`, `128` and `64` for binary, ternary and quaternary mixtures (as in [1]), and `32` for more components (smaller for more than six components, at most ~5⋅10⁵ grid points). Ignored if `diagram` is given.
- `equilibrium`: `:auto` considers liquid and vapour phases, `:lle` only liquid phases. Ignored if `diagram` is given.
- `polish`: `TPFlashMethod` used to refine the CEM phases (initialized with the CEM result), or `nothing` to return the CEM result. The refined result is only accepted if its phases are in equilibrium and it has a lower Gibbs energy than the CEM result. If `MichelsenTPFlash` is given and the CEM finds more than two phases, a `MultiPhaseTPFlash` with the same tolerances is used. Results with solid phases are not refined.
- `threaded`: evaluate the grid in parallel.
- `diagram`: a precomputed `CEMDiagram`.

Pure solid phases are considered for `CompositeModel`s with a `SolidHfus` solid model. Solid phases are returned with a `NaN` volume.

## References
1. Göttl, Q., Rosen, N., & Burger, J. (2025). Convex envelope method for T, p flash calculations for mixtures with an arbitrary number of components and arbitrary aggregate states. Computers & Chemical Engineering, 203, 109326. [doi:10.1016/j.compchemeng.2025.109326](https://doi.org/10.1016/j.compchemeng.2025.109326)
"""
Base.@kwdef struct CEMFlash <: TPFlashMethod
    δ::Union{Int,Nothing} = nothing
    equilibrium::Symbol = :auto
    polish::Union{Nothing,TPFlashMethod} = MichelsenTPFlash(second_order = true)
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
    ws = Solvers.EnvelopeWorkspace(env.lp)
    zf = collect(Float64,z)
    Solvers.envelope_reset!(ws,env.lp,zf)
    Solvers.envelope_lp!(ws,env.lp,zf)
    result = cem_flash_result(env,fluid,p,T,zf,ws,stable)
    return cem_polish(fluidmodel,fluid,p,T,zf,result,method.polish)
end

cem_polish(model,fluid,p,T,z,result,::Nothing) = result

function cem_polish(model,fluid,p,T,z,result,method::TPFlashMethod)
    #single phase, or solid phases present (the local methods do not handle solids)
    (numphases(result) < 2 || !all(isfinite,result.volumes)) && return result
    m = cem_polish_method(method,result)
    m === nothing && return result
    polished = try
        tp_flash_impl(model,p,T,z,m)
    catch
        return result
    end
    #the discrete envelope lies above the exact one: a correct refinement can only lower g
    g0,g1 = result.data.g,cem_equilibrium_gibbs(fluid,p,T,z,polished)
    return g1 <= g0 + 1e-10*(1 + abs(g0)) ? polished : result
end

#=
reduced gibbs energy of a flash result, with the same reference as the envelope, NaN if its phases are not in equilibrium.
the local methods may return without converging, and their gibbs energies do not always use the same reference.
=#
function cem_equilibrium_gibbs(fluid,p,T,z,res)
    cache = ∂lnϕ_cache(fluid,p,T,z,Val{false}())
    g = 0.0
    lnf1 = nothing
    for j in 1:numphases(res)
        β,x = res.fractions[j],res.compositions[j]
        iszero(β) && continue
        phase = j == res.data.vapour_idx ? :vapour : :liquid
        gj,v = modified_gibbs(fluid,p,T,x,phase,NaN)
        g += β*gj
        lnϕ,_ = modified_lnϕ(fluid,p,T,x,cache;phase,vol0 = v)
        if lnf1 === nothing
            lnf1 = lnϕ .+ log.(x)
        elseif !all(i -> iszero(z[i]) || abs(lnϕ[i] + log(x[i]) - lnf1[i]) <= CEM_LNF_TOL,eachindex(z))
            return NaN
        end
    end
    return g
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
