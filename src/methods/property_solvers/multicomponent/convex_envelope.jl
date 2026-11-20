#=
Convex envelope method (CEM)

Göttl, Rosen, Burger (2025), Comput. Chem. Eng. 203, 109326. doi:10.1016/j.compchemeng.2025.109326

The reduced Gibbs energy of every aggregate state is evaluated on a lattice of the composition simplex
(step 1/δ). The phase split of a feed z is the facet of the lower convex envelope of those points above z,
found with the linear program in Solvers.LowerEnvelope.

Points of the optimal basis belong to the same phase when they have the same aggregate state (fluid points
with almost equal volumes are also merged: near-critical region) and are lattice neighbours, or are connected
through the stable single phase region (see cem_samephase).
=#

@enum CEMState::Int8 CEM_LIQUID = 1 CEM_VAPOUR CEM_SOLID

cem_state_symbol(s) = s == CEM_LIQUID ? :liquid : (s == CEM_VAPOUR ? :vapour : :solid)

#N > 4: at most 32, and at most ~5e5 grid points (the grid grows as δ^(N-1))
cem_default_δ(N) = N == 2 ? 256 : (N == 3 ? 128 : (N == 4 ? 64 : max(1,findlast(δ -> binomial(δ + N - 1,N - 1) <= 500_000,1:32))))

#grid columns: one per lattice point and aggregate state, in the (tree) order of lp
struct CEMEnvelope
    δ::Int
    lattice::Matrix{Int32} #N × npoints, lattice coordinates (x = k/δ)
    lp::Solvers.LowerEnvelope #g: reduced gibbs energy, minus the linear function x⋅gpure
    gpure::Vector{Float64}
    point::Vector{Int32} #lattice point of each column
    state::Vector{CEMState}
    volume::Vector{Float64}
    point_column::Vector{Int32} #column with the lowest g at each lattice point
    point_ptr::Vector{Int32} #columns of lattice point k: point_cols[point_ptr[k]:point_ptr[k+1]-1]
    point_cols::Vector{Int32}
    λ_tol::Float64 #phase amount tolerance
    vol_rtol::Float64 #fluid points with this relative volume difference are the same phase
end

Base.length(env::CEMEnvelope) = size(env.lattice,1)

#model used to evaluate the fluid phases, solid model (or nothing)
function cem_models(model,p,T,z,equilibrium)
    solid = nothing
    fluid = model
    if model isa CompositeModel && model.solid !== nothing
        solid = model.solid
        fluid = model.fluid
        if !(solid isa SolidHfusModel)
            throw(ArgumentError("CEM: only `SolidHfus` solid models are supported (pure solid phases), got $(typeof(solid))"))
        end
    end
    cached_fluid = __tpflash_cache_model(fluid,p,T,z,is_lle(equilibrium) ? :lle : :vle) #PTFlashWrapper does not know :auto
    return fluid,cached_fluid,solid
end

#Δμᵢ/RT (pure solid - pure liquid) of each component, NaN if the component does not form a solid
function cem_solid_Δμ(model,solid,p,T)
    Δμ = fill(NaN,length(model))
    solid === nothing && return Δμ
    μ = chemical_potential(solid,p,T,ones(length(solid)))
    idx = _mapping_split(model)
    ν = _mapping_fractions(model)
    for j in eachindex(idx)
        length(idx[j]) == 1 || throw(ArgumentError("CEM: only pure solids are supported, solid $(component_list(solid)[j]) is formed by several fluid components"))
        Δμ[only(idx[j])] = μ[j]/(only(ν[j])*Rgas(solid)*T)
    end
    return Δμ
end

#returns (gl,vl,gv,vv), NaN g if the column does not exist (lle, single EoS root of the other state, missing root)
function cem_eval_point(fluid,p,T,x,lle)
    if fluid isa PTFlashWrapper
        gl,vl = modified_gibbs(fluid,p,T,x,:liquid,NaN)
        gv,vv = lle ? (NaN,NaN) : modified_gibbs(fluid,p,T,x,:vapour,NaN)
    else
        vl = volume(fluid,p,T,x,phase = :liquid)
        gl,_ = modified_gibbs(fluid,p,T,x,:liquid,vl)
        gv,vv = NaN,NaN
        if !lle
            vv = volume(fluid,p,T,x,phase = :vapour)
            if isapprox(vl,vv,rtol = 1e-8)
                #only one root: store it with the state it actually has
                gv,vv = NaN,NaN
                is_vapour(identify_phase(fluid,p,T,x,vol = vl)) && ((gl,vl,gv,vv) = (gv,vv,gl,vl))
            else
                gv,_ = modified_gibbs(fluid,p,T,x,:vapour,vv)
            end
        end
    end
    #a missing root (e.g. no liquid root above the critical point) gives v = Inf and a meaningless finite g
    return (isfinite(vl) ? gl : NaN),vl,(isfinite(vv) ? gv : NaN),vv
end

function cem_chunks(n,threaded)
    #more chunks than threads: the cost per lattice point is very uneven (load balancing)
    nt = threaded ? 8*Threads.nthreads() : 1
    nt = max(1,min(nt,n ÷ 256))
    #contiguous blocks: consecutive lattice points are neighbours (warm starts)
    bounds = round.(Int,range(0,n,length = nt + 1))
    return [(bounds[i] + 1):bounds[i + 1] for i in 1:nt]
end

function cem_foreach_chunk(f,n,threaded)
    chunks = cem_chunks(n,threaded)
    if length(chunks) == 1
        return [f(chunks[1])]
    else
        tasks = [StableTasks.@spawn f(chunk) for chunk in chunks]
        return fetch.(tasks)
    end
end

function cem_eval_chunk!(gl,vl,gv,vv,fluid,p,T,K,δ,lle,idx)
    x = zeros(size(K,1))
    for k in idx
        x .= @view(K[:,k]) ./ δ
        gl[k],vl[k],gv[k],vv[k] = cem_eval_point(fluid,p,T,x,lle)
    end
    return nothing
end

#evaluates the Gibbs energy of all aggregate states on the composition lattice and builds the envelope
function CEMEnvelope(model,p,T,method)
    N = length(model)
    δ = something(method.δ,cem_default_δ(N))
    equilibrium,threaded = method.equilibrium,method.threaded
    δ < 1 && throw(DomainError(δ,"CEM: δ must be a positive integer"))
    _,fluid,solid = cem_models(model,p,T,fill(1/N,N),equilibrium)
    lle = is_lle(equilibrium)
    K = Solvers.simplex_lattice(N,δ)
    M = size(K,2)
    gl,vl,gv,vv = zeros(M),zeros(M),zeros(M),zeros(M)
    cem_foreach_chunk(idx -> cem_eval_chunk!(gl,vl,gv,vv,fluid,p,T,K,δ,lle,idx),M,threaded)
    #compact into columns
    P,G,S,V = Int32[],Float64[],CEMState[],Float64[]
    function addcol!(k,g,s,v)
        isfinite(g) || return nothing
        push!(P,k);push!(G,g);push!(S,s);push!(V,v)
        return nothing
    end
    for k in 1:M
        addcol!(k,gl[k],CEM_LIQUID,vl[k])
        addcol!(k,gv[k],CEM_VAPOUR,vv[k])
    end
    pure = [Solvers.simplex_lattice_rank(δ .* (1:N .== i),δ) for i in 1:N]
    for (k,Δμ) in zip(pure,cem_solid_Δμ(model,solid,p,T))
        #pure solid, relative to the pure liquid (or the only fluid root)
        addcol!(k,(isfinite(gl[k]) ? gl[k] : gv[k]) + Δμ,CEM_SOLID,NaN)
    end
    X = @view(K[:,P]) ./ δ
    #remove the linear part (affine transformations do not change the envelope)
    gpure = [minimum((G[c] for c in eachindex(G) if P[c] == k),init = Inf) for k in pure]
    all(isfinite,gpure) || error("CEM: Gibbs energy evaluation failed for a pure component at p = $p, T = $T")
    G .-= X'*gpure
    lp,perm = Solvers.lower_envelope(X,G,Int8.(S))
    point = P[perm]
    point_column = zeros(Int32,M)
    for c in eachindex(point)
        c0 = point_column[point[c]]
        (c0 == 0 || lp.g[c] < lp.g[c0]) && (point_column[point[c]] = c)
    end
    point_cols = Int32.(sortperm(point))
    point_ptr = zeros(Int32,M + 1)
    point_ptr[1] = 1
    for k in point
        point_ptr[k + 1] += 1
    end
    cumsum!(point_ptr,point_ptr)
    return CEMEnvelope(δ,K,lp,gpure,point,S[perm],V[perm],point_column,point_ptr,point_cols,method.λ_tol,method.vol_rtol)
end

function cem_add_point_columns!(cand,env::CEMEnvelope,k)
    for i in env.point_ptr[k]:(env.point_ptr[k + 1] - 1)
        push!(cand,env.point_cols[i])
    end
    return cand
end

#columns at lattice point k and its neighbours
function cem_add_neighbourhood!(cand,env::CEMEnvelope,k,kbuf)
    K = env.lattice
    N = size(K,1)
    cem_add_point_columns!(cand,env,k)
    for j in 1:N
        K[j,k] == 0 && continue
        for i in 1:N
            i == j && continue
            kbuf .= @view(K[:,k])
            kbuf[j] -= 1
            kbuf[i] += 1
            cem_add_point_columns!(cand,env,Solvers.simplex_lattice_rank(kbuf,env.δ))
        end
    end
    return cand
end

#candidates for the local LP of lattice point k: neighbourhood of k and of the current basis
function cem_local_candidates!(cand,kbuf,ws::Solvers.EnvelopeWorkspace,env::CEMEnvelope,k)
    empty!(cand)
    cem_add_neighbourhood!(cand,env,k,kbuf)
    for c in ws.B
        cem_add_neighbourhood!(cand,env,env.point[c],kbuf)
    end
    return cand
end

## phases

#columns a and b can belong to the same phase (aggregate state)
function cem_compatible(env::CEMEnvelope,a,b)
    sa,sb = env.state[a],env.state[b]
    (sa == CEM_SOLID || sb == CEM_SOLID) && return a == b
    sa == sb && return true
    return isapprox(env.volume[a],env.volume[b],rtol = env.vol_rtol) #near critical: the roots merge
end

#stability (single phase) of lattice points, from a diagram or lazily from the LP of each point
struct CEMStability
    env::CEMEnvelope
    point_facet::Union{Nothing,Vector{Int32}}
    ws::Solvers.EnvelopeWorkspace
    memo::Dict{Int,Bool}
end

CEMStability(env::CEMEnvelope,point_facet = nothing) = CEMStability(env,point_facet,Solvers.EnvelopeWorkspace(env.lp),Dict{Int,Bool}())

function (stable::CEMStability)(k)
    stable.point_facet !== nothing && return iszero(stable.point_facet[k])
    return get!(stable.memo,k) do
        env,ws = stable.env,stable.ws
        z = env.lp.x[:,env.point_column[k]]
        Solvers.envelope_reset!(ws,env.lp,z)
        Solvers.envelope_lp!(ws,env.lp,z)
        numphases(env,ws) == 1
    end
end

#=
Two hull vertices of the same aggregate state are the same phase if they are lattice neighbours [1].
At grid resolution, a binodal that runs almost parallel to a lattice direction (e.g. close to an edge) has hull vertices that are
several lattice steps apart. They are also the same phase if the segment between them stays in the stable single phase region,
i.e. a stable lattice point of the same state lies on or next to every lattice step of the segment. For a real phase split, the
segment crosses the miscibility gap.
=#
function cem_samephase(env::CEMEnvelope,a,b,stable = nothing)
    cem_compatible(env,a,b) || return false
    K = env.lattice
    pa,pb = env.point[a],env.point[b]
    dist,n = 0,0
    for i in 1:size(K,1)
        dk = abs(K[i,pa] - K[i,pb])
        dist += dk
        n = max(n,dk)
    end
    dist <= 2 && return true
    (stable === nothing || n > max(4,env.δ ÷ 8)) && return false
    N = size(K,1)
    x = zeros(N)
    kq = zeros(Int32,N)
    kn = zeros(Int32,N)
    for j in 1:(n - 1)
        for i in 1:N
            x[i] = K[i,pa] + (j/n)*(K[i,pb] - K[i,pa])
        end
        Solvers.simplex_lattice_round!(kq,x,env.δ)
        ok = cem_stable_compatible(env,Solvers.simplex_lattice_rank(kq,env.δ),a,stable)
        for jj in 1:N, ii in 1:N
            (ok || ii == jj || kq[jj] == 0) && continue
            kn .= kq
            kn[jj] -= 1
            kn[ii] += 1
            ok = cem_stable_compatible(env,Solvers.simplex_lattice_rank(kn,env.δ),a,stable)
        end
        ok || return false
    end
    return true
end

cem_stable_compatible(env,q,a,stable) = stable(q) && cem_compatible(env,env.point_column[q],a)

#cluster index of each column in cols, number of clusters
function cem_clusters(env::CEMEnvelope,cols,stable = nothing)
    n = length(cols)
    label = collect(1:n) #union-find parent, always a smaller index
    root(i) = (while label[i] != i; i = label[i]; end; i)
    for a in 1:n, b in (a + 1):n
        ra,rb = root(a),root(b)
        ra != rb && cem_samephase(env,cols[a],cols[b],stable) && (label[max(ra,rb)] = min(ra,rb))
    end
    #parents come first: number the roots, children take the label of their parent
    np = 0
    for i in 1:n
        label[i] = label[i] == i ? (np += 1) : label[label[i]]
    end
    return label,np
end

function numphases(env::CEMEnvelope,ws::Solvers.EnvelopeWorkspace,stable = nothing)
    count(>(env.λ_tol),ws.λ) <= 1 && return 1
    return last(cem_clusters(env,ws.B[ws.λ .> env.λ_tol],stable))
end

#phases (sorted by volume, solids last) of the current basis with weights λ
function cem_phases(env::CEMEnvelope,fluid,p,T,B,λ,stable = nothing)
    support = findall(>(env.λ_tol),λ)
    cols = B[support]
    w = λ[support]
    labels,np = cem_clusters(env,cols,stable)
    N = length(env)
    comps = [zeros(N) for _ in 1:np]
    β = zeros(np)
    states = fill(CEM_LIQUID,np)
    wmax = zeros(np)
    for (j,c) in enumerate(cols)
        l = labels[j]
        β[l] += w[j]
        comps[l] .+= w[j] .* @view(env.lp.x[:,c])
        if w[j] > wmax[l]
            wmax[l] = w[j]
            states[l] = env.state[c]
        end
    end
    vols = zeros(np)
    for l in 1:np
        comps[l] ./= β[l]
        s = states[l]
        vols[l] = s == CEM_SOLID ? NaN : volume(fluid,p,T,comps[l],phase = cem_state_symbol(s))
    end
    idx = sortperm(vols) #NaN (solids) last
    return comps[idx],β[idx],vols[idx],states[idx]
end

function FlashResult(env::CEMEnvelope,fluid,p,T,z,ws::Solvers.EnvelopeWorkspace,stable = nothing)
    comps,β,vols,states = cem_phases(env,fluid,p,T,ws.B,ws.λ,stable)
    g = dot(ws.λ,@view(env.lp.g[ws.B])) + dot(z,env.gpure)
    vapour_idx = something(findlast(==(CEM_VAPOUR),states),-1) #largest volume
    return FlashResult(comps,β,vols,FlashData(p,T,g,vapour_idx))
end

## phase diagram

"""
    CEMPhaseDiagram(;δ = nothing, equilibrium = :auto, polish = nothing, threaded = true, λ_tol = 1e-12, vol_rtol = 0.1, lnf_tol = 1e-6)

Method to calculate the complete phase diagram at constant pressure and temperature with the convex envelope method (CEM) [1].

The reduced Gibbs energy of every aggregate state (liquid, vapour and, for `CompositeModel` with a `SolidHfus` solid model, pure solids) is evaluated on a lattice of the composition simplex with step `1/δ`.
The lower convex envelope of those points is then determined for every lattice point. Its facets that connect points of different phases are the multiphase regions (tie-lines, tie-triangles, …).
Instead of building the complete convex hull, the envelope is obtained by solving a linear program for every lattice point, warm-started from the previous point. The tangent plane distances needed by the simplex method are evaluated with a tree that discards most of the grid with a rigorous bound.

## Keyword arguments
- `δ`: number of divisions of the composition axes. By default, `256`, `128` and `64` for binary, ternary and quaternary mixtures (as in [1]), and `32` for more components (smaller for more than six components, at most ~5⋅10⁵ grid points). The number of grid points is `binomial(δ + N - 1, N - 1)`.
- `equilibrium`: `:auto` considers liquid and vapour phases, `:lle` only liquid phases.
- `polish`: `nothing` (default) or a `TPFlashMethod` (e.g. `MichelsenTPFlash()` or `MultiPhaseTPFlash()`). If given, the phases of every multiphase facet are refined with that method, starting from the CEM phases. The refined result is only accepted if its phases are in equilibrium and it lowers the Gibbs energy. Facets with solid phases are not refined. Close to critical points, `MichelsenTPFlash(second_order = true)` converges more reliably.
- `threaded`: evaluate the grid and the envelope in parallel.
- `λ_tol`: vertices of a facet with a smaller weight are not part of the phase split.
- `vol_rtol`: liquid and vapour grid points whose volumes differ by less than this relative tolerance can belong to the same phase (near critical region).
- `lnf_tol`: tolerance of the log fugacities for a polished result to be accepted.

## References
1. Göttl, Q., Rosen, N., & Burger, J. (2025). Convex envelope method for T, p flash calculations for mixtures with an arbitrary number of components and arbitrary aggregate states. Computers & Chemical Engineering, 203, 109326. [doi:10.1016/j.compchemeng.2025.109326](https://doi.org/10.1016/j.compchemeng.2025.109326)
"""
Base.@kwdef struct CEMPhaseDiagram <: ThermodynamicMethod
    δ::Union{Int,Nothing} = nothing
    equilibrium::Symbol = :auto
    polish::Union{Nothing,TPFlashMethod} = nothing
    threaded::Bool = true
    λ_tol::Float64 = 1e-12
    vol_rtol::Float64 = 0.1
    lnf_tol::Float64 = 1e-6
end

#phases (sorted by volume, solids last with NaN volume). Without polishing, each phase composition is the mean of its facet vertices
struct CEMFacet
    columns::Vector{Int32} #grid columns of the facet vertices
    μ::Vector{Float64} #tangent plane μᵢ/RT, in the reference state of the model
    compositions::Vector{Vector{Float64}}
    volumes::Vector{Float64}
    states::Vector{Symbol}
end

struct CEMDiagram
    components::Vector{String}
    p::Float64
    T::Float64
    equilibrium::Symbol
    envelope::CEMEnvelope
    facets::Vector{CEMFacet} #multiphase facets
    point_facet::Vector{Int32} #facet of each grid point, 0: stable single phase
end

Base.show(io::IO,env::CEMEnvelope) = print(io,"CEMEnvelope(δ = ",env.δ,", ",size(env.lattice,2)," grid points, ",length(env.lp.g)," columns)")
Base.show(io::IO,d::CEMDiagram) = print(io,"CEMDiagram(",length(d.components)," components, p = ",d.p," Pa, T = ",d.T," K, δ = ",d.envelope.δ,", ",length(d.facets)," multiphase facets)")

function Base.show(io::IO,::MIME"text/plain",d::CEMDiagram)
    env = d.envelope
    nps = [length(f.compositions) for f in d.facets]
    print(io,"CEMDiagram(",length(d.components)," components, p = ",d.p," Pa, T = ",d.T," K, δ = ",env.δ,", ",size(env.lattice,2)," grid points)")
    print(io,"\n multiphase facets: ",length(d.facets))
    for np in sort(unique(nps))
        print(io,"\n  ",np,"-phase: ",count(isequal(np),nps))
    end
    print(io,"\n single phase grid points: ",count(iszero,d.point_facet))
end

"""
    phase_diagram(model, p, T, method = CEMPhaseDiagram())
    phase_diagram(model, p, T; kwargs...)

Calculates the complete phase diagram of `model` at pressure `p` `[Pa]` and temperature `T` `[K]`, over the whole composition space.
Returns a `CEMDiagram`, containing all multiphase regions (as facets of the convex envelope of the Gibbs energy) and the phase classification of every grid point.

Currently, the only method available is [`CEMPhaseDiagram`](@ref). The keyword arguments are passed to it.

## Example
```julia
model = UNIFAC(["water","hexane","ethanol"])
diagram = phase_diagram(model, 1e5, 298.15, CEMPhaseDiagram(δ = 64))
diagram.facets[1].compositions #compositions of the phases of the first multiphase facet
tp_flash2(model, 1e5, 298.15, [0.4,0.4,0.2], CEMFlash(diagram = diagram)) #flash with the stored envelope
```
"""
function phase_diagram(model::EoSModel,p,T;kwargs...)
    return phase_diagram(model,p,T,CEMPhaseDiagram(;kwargs...))
end

function phase_diagram(model::EoSModel,_p,_T,method::CEMPhaseDiagram)
    p,T = ustrip(_p,pressure),ustrip(_T,temperature)
    p,T = float(primalval(p)),float(primalval(T))
    N = length(model)
    N < 2 && throw(ArgumentError("phase_diagram: at least two components are required"))
    env = CEMEnvelope(model,p,T,method)
    fluidmodel,fluid,_ = cem_models(model,p,T,fill(1/N,N),method.equilibrium)
    M = size(env.lattice,2)
    walk_facet = zeros(Int32,M)
    #each chunk walks its lattice points in order, warm-starting every LP from the previous solution
    results = cem_foreach_chunk(idx -> cem_walk!(walk_facet,env,idx),M,method.threaded)
    #merge the facets found by every chunk
    ids = Dict{Vector{Int32},Int}()
    walk_cols = Vector{Int32}[]
    for (idx,list) in results
        newid = map(list) do key
            get!(ids,key) do
                push!(walk_cols,key)
                length(walk_cols)
            end
        end
        for k in idx
            walk_facet[k] != 0 && (walk_facet[k] = newid[walk_facet[k]])
        end
    end
    facet_cols,point_facet = cem_merge_phases(env,walk_cols,walk_facet)
    stable = CEMStability(env,point_facet)
    ws = Solvers.EnvelopeWorkspace(env.lp)
    w = fill(1/N,N)
    facets = map(facet_cols) do cols
        ws.B .= cols
        z = env.lp.x[:,cols]*w
        Solvers.envelope_refactor!(ws,env.lp,z)
        μ = Solvers.envelope_tangent_plane(ws,env.lp) .+ env.gpure
        comps,_,vols,states = cem_phases(env,fluid,p,T,cols,w,stable)
        f = CEMFacet(cols,μ,comps,vols,cem_state_symbol.(states))
        method.polish === nothing ? f : cem_polish_facet(fluidmodel,p,T,z,f,env,fluid,ws,stable,method.polish,method.lnf_tol)
    end
    return CEMDiagram(component_list(model),p,T,method.equilibrium,env,facets,point_facet)
end

conodes(d::CEMDiagram) = [f.compositions for f in d.facets]

#second pass over the multiphase facets found by the walk: vertices that are the same phase but not lattice neighbours
#(see cem_samephase) are merged, points and facets that become single phase are removed.
function cem_merge_phases(env::CEMEnvelope,facet_cols,point_facet)
    N = length(env)
    stable = CEMStability(env,copy(point_facet))
    labels = [first(cem_clusters(env,cols,stable)) for cols in facet_cols]
    Binv = [inv(env.lp.x[:,cols]) for cols in facet_cols]
    λ = zeros(N)
    newpf = zeros(Int32,length(point_facet))
    used = falses(length(facet_cols))
    for k in eachindex(point_facet)
        f = point_facet[k]
        f == 0 && continue
        mul!(λ,Binv[f],@view(env.lp.x[:,env.point_column[k]]))
        l = labels[f]
        i0 = findfirst(>(env.λ_tol),λ)
        if any(i -> λ[i] > env.λ_tol && l[i] != l[i0],1:N)
            newpf[k] = f
            used[f] = true
        end
    end
    newid = cumsum(used)
    for k in eachindex(newpf)
        newpf[k] != 0 && (newpf[k] = newid[newpf[k]])
    end
    return facet_cols[used],newpf
end

#solves the LP for the lattice points idx in order, warm-starting from the previous solution.
#optimality (all reduced costs ≥ 0) only depends on the basis, so a proven basis resolves without any pricing:
# - the next point, if it lies inside the basis simplex (λ = B⁻¹z ≥ 0),
# - all vertices of the basis: their columns touch the envelope, so they are stable single phase points.
function cem_walk!(point_facet,env::CEMEnvelope,idx)
    ws = Solvers.EnvelopeWorkspace(env.lp)
    keys = Dict{Vector{Int32},Int}()
    list = Vector{Int32}[]
    proven = Set{Vector{Int32}}()
    resolved = falses(size(env.lattice,2))
    key = similar(ws.B)
    cand,kbuf = Int32[],similar(ws.B)
    z = zeros(length(env))
    λ = ws.λ
    current_proven = false
    for k in idx
        resolved[k] && continue
        c = env.point_column[k]
        z .= @view(env.lp.x[:,c])
        inside = false
        if current_proven
            mul!(λ,ws.Binv,z)
            inside = all(>=(-1e-12),λ)
            inside && (λ .= max.(λ,0.0))
        end
        if !inside
            if !current_proven
                Solvers.envelope_reset!(ws,env.lp,z)
            else
                cem_local_candidates!(cand,kbuf,ws,env,k) #before the crash: neighbourhood of the previous basis
                Solvers.envelope_crash!(ws,env.lp,c,z)
                Solvers.envelope_lp!(ws,env.lp,z,cand)
            end
            sort!(key .= ws.B)
            if !(key in proven)
                Solvers.envelope_lp!(ws,env.lp,z)
                sort!(key .= ws.B)
                push!(proven,copy(key))
            end
            for cB in ws.B
                resolved[env.point[cB]] = true
            end
            current_proven = true
        end
        if numphases(env,ws) > 1
            sort!(key .= ws.B)
            id = get(keys,key,0)
            if id == 0
                fkey = copy(key) #key is a reused buffer
                push!(list,fkey)
                id = keys[fkey] = length(list)
            end
            point_facet[k] = id
        end
    end
    return idx,list
end

function cem_polish_facet(model,p,T,z,f::CEMFacet,env,fluid,ws,stable,method,lnf_tol)
    ws.λ .= 1/length(z)
    res = FlashResult(env,fluid,p,T,z,ws,stable)
    res2 = cem_polish(model,fluid,p,T,z,res,method,lnf_tol)
    res2 === res && return f
    np = numphases(res2,true)
    np == length(f.compositions) || return f
    idx = [i for i in 1:numphases(res2) if is_active_phase(res2,i)]
    vapour_idx = res2.data.vapour_idx
    states = [i == vapour_idx ? :vapour : :liquid for i in idx]
    return CEMFacet(f.columns,f.μ,[collect(res2.compositions[i]) for i in idx],res2.volumes[idx],states)
end

export phase_diagram, CEMPhaseDiagram
