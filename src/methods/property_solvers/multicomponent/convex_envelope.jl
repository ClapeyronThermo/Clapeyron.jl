#=
Convex envelope method (CEM)

Göttl, Rosen, Burger (2025), Comput. Chem. Eng. 203, 109326. doi:10.1016/j.compchemeng.2025.109326

The reduced Gibbs energy of every aggregate state is evaluated on a lattice of the composition simplex
(step 1/δ). The phase split of a feed z is the facet of the lower convex envelope of those points above z:

    min Σₖ λₖ gₖ   s.t.  Σₖ λₖ xₖ = z,  λ ≥ 0

This linear program is solved with a revised simplex method on an N×N basis instead of building the hull.
Pricing (the tangent plane distance of every grid point to the current tangent plane) is accelerated
with a kd-tree that stores, for every node, an affine under-estimator of g. This gives a rigorous lower
bound of the reduced cost on the node, so most of the grid is never touched.

Points of the optimal basis belong to the same phase when they have the same aggregate state (fluid points
with almost equal volumes are also merged: near-critical region) and are lattice neighbours, or are connected
through the stable single phase region (see cem_samephase).
=#

const CEM_LIQUID = Int8(1)
const CEM_VAPOUR = Int8(2)
const CEM_SOLID = Int8(3)
const CEM_RTOL = 1e-10 #reduced cost tolerance
const CEM_λTOL = 1e-12 #phase amount tolerance
const CEM_VOL_RTOL = 0.1 #fluid points with this relative volume difference are the same phase

cem_state_symbol(s) = s == CEM_LIQUID ? :liquid : (s == CEM_VAPOUR ? :vapour : :solid)

cem_default_δ(N) = N == 2 ? 256 : (N == 3 ? 128 : (N == 4 ? 64 : 32))

#grid + kd tree over the grid columns. columns are stored in tree order.
struct CEMEnvelope
    δ::Int
    lattice::Matrix{Int32} #N × npoints, lattice coordinates (x = k/δ)
    x::Matrix{Float64} #N × ncolumns
    g::Vector{Float64} #reduced gibbs energy, minus the linear function x⋅gpure
    gpure::Vector{Float64}
    point::Vector{Int32} #lattice point of each column
    state::Vector{Int8}
    volume::Vector{Float64}
    point_column::Vector{Int32} #column with the lowest g at each lattice point
    point_ptr::Vector{Int32} #columns of lattice point k: point_cols[point_ptr[k]:point_ptr[k+1]-1]
    point_cols::Vector{Int32}
    vertex::Vector{Int32} #point_column of the pure components
    node_range::Matrix{Int32} #2 × nnodes
    node_child::Matrix{Int32} #2 × nnodes, 0 for leaves
    node_lo::Matrix{Float64} #(N-1) × nnodes, bounding box in reduced coordinates (x₁ … xₙ₋₁)
    node_hi::Matrix{Float64}
    node_b::Matrix{Float64} #slope of the affine under-estimator
    node_e::Vector{Float64} #minimum of g - b⋅x on the node
end

Base.length(env::CEMEnvelope) = size(env.x,1)

function cem_lattice(N,δ)
    K = Matrix{Int32}(undef,N,binomial(δ + N - 1,N - 1))
    k = zeros(Int32,N)
    _cem_lattice!(K,k,1,δ,Ref(0))
    return K
end

function _cem_lattice!(K,k,i,rem,j)
    if i == length(k)
        k[i] = rem
        j[] += 1
        K[:,j[]] .= k
        return nothing
    end
    for ki in rem:-1:0
        k[i] = ki
        _cem_lattice!(K,k,i + 1,rem - ki,j)
    end
    return nothing
end

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
    N = length(model)
    Δμ = fill(NaN,N)
    solid === nothing && return Δμ
    solid_comps = component_list(solid)
    μ = chemical_potential(solid,p,T,ones(length(solid)))
    fluid_comps = component_list(model)
    for i in 1:N
        j = findfirst(isequal(fluid_comps[i]),solid_comps)
        j === nothing || (Δμ[i] = μ[j]/(Rgas(solid)*T))
    end
    return Δμ
end

#returns (gl,vl,gv,vv), NaN g if the column does not exist (lle, or single EoS root of the other state)
function cem_eval_point(fluid,p,T,x,lle)
    nan = NaN
    if fluid isa PTFlashWrapper
        gl,vl = modified_gibbs(fluid,p,T,x,:liquid,nan)
        lle && return gl,vl,nan,nan
        gv,vv = modified_gibbs(fluid,p,T,x,:vapour,nan)
        return gl,vl,gv,vv
    end
    vl = volume(fluid,p,T,x,phase = :liquid)
    gl,_ = modified_gibbs(fluid,p,T,x,:liquid,vl)
    lle && return gl,vl,nan,nan
    vv = volume(fluid,p,T,x,phase = :vapour)
    if isapprox(vl,vv,rtol = 1e-8)
        #only one root: store it with the state it actually has
        if is_vapour(identify_phase(fluid,p,T,x,vol = vl))
            return nan,nan,gl,vl
        else
            return gl,vl,nan,nan
        end
    end
    gv,_ = modified_gibbs(fluid,p,T,x,:vapour,vv)
    return gl,vl,gv,vv
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

"""
    cem_envelope(model,p,T,δ,equilibrium = :auto,threaded = true)

Evaluates the Gibbs energy of all aggregate states on the composition lattice and builds the search tree.
"""
function cem_envelope(model,p,T,δ,equilibrium = :auto,threaded = true)
    N = length(model)
    δ < 1 && throw(DomainError(δ,"CEM: δ must be a positive integer"))
    z0 = fill(1/N,N)
    _,fluid,solid = cem_models(model,p,T,z0,equilibrium)
    lle = is_lle(equilibrium)
    K = cem_lattice(N,δ)
    M = size(K,2)
    gl,vl,gv,vv = zeros(M),zeros(M),zeros(M),zeros(M)
    cem_foreach_chunk(idx -> cem_eval_chunk!(gl,vl,gv,vv,fluid,p,T,K,δ,lle,idx),M,threaded)
    Δμs = cem_solid_Δμ(model,solid,p,T)
    #compact into columns
    Xv = Float64[]
    G,P,S,V = Float64[],Int32[],Int8[],Float64[]
    function addcol!(k,g,s,v)
        isfinite(g) || return nothing
        append!(Xv,@view(K[:,k]))
        push!(G,g);push!(P,k);push!(S,s);push!(V,v)
        return nothing
    end
    for k in 1:M
        addcol!(k,gl[k],CEM_LIQUID,vl[k])
        addcol!(k,gv[k],CEM_VAPOUR,vv[k])
        i = findfirst(isequal(δ),@view(K[:,k]))
        if i !== nothing && isfinite(Δμs[i])
            #pure solid, relative to the pure liquid (or the only fluid root)
            addcol!(k,(isfinite(gl[k]) ? gl[k] : gv[k]) + Δμs[i],CEM_SOLID,NaN)
        end
    end
    X = reshape(Xv,N,length(G)) ./ δ
    #remove the linear part (affine transformations do not change the envelope)
    gpure = fill(Inf,N)
    for c in eachindex(G)
        i = findfirst(isequal(δ),@view(K[:,P[c]]))
        i !== nothing && (gpure[i] = min(gpure[i],G[c]))
    end
    all(isfinite,gpure) || throw(error("CEM: Gibbs energy evaluation failed for a pure component at p = $p, T = $T"))
    for c in eachindex(G)
        G[c] -= dot(@view(X[:,c]),gpure)
    end
    return cem_build_envelope(δ,K,X,G,gpure,P,S,V)
end

## kd tree
struct CEMTreeBuilder
    x::Matrix{Float64}
    g::Vector{Float64}
    state::Vector{Int8}
    perm::Vector{Int32}
    range::Vector{Int32}
    child::Vector{Int32}
    lo::Vector{Float64}
    hi::Vector{Float64}
    b::Vector{Float64}
    e::Vector{Float64}
    leafsize::Int
end

function cem_build_envelope(δ,K,X,G,gpure,P,S,V)
    N,nc = size(X)
    tb = CEMTreeBuilder(X,G,S,collect(Int32,1:nc),Int32[],Int32[],Float64[],Float64[],Float64[],Float64[],16)
    cem_build_node!(tb,1,nc)
    perm = tb.perm
    nn = length(tb.e)
    x = X[:,perm]
    g = G[perm]
    point = P[perm]
    state = S[perm]
    vol = V[perm]
    M = size(K,2)
    point_column = zeros(Int32,M)
    for c in 1:nc
        k = point[c]
        c0 = point_column[k]
        if c0 == 0 || g[c] < g[c0]
            point_column[k] = c
        end
    end
    vertex = zeros(Int32,N)
    for k in 1:M
        i = findfirst(isequal(δ),@view(K[:,k]))
        i !== nothing && (vertex[i] = point_column[k])
    end
    point_cols = Int32.(sortperm(point))
    point_ptr = zeros(Int32,M + 1)
    point_ptr[1] = 1
    for c in 1:nc
        point_ptr[point[c] + 1] += 1
    end
    cumsum!(point_ptr,point_ptr)
    nr = max(N - 1,1)
    return CEMEnvelope(δ,K,x,g,gpure,point,state,vol,point_column,point_ptr,point_cols,vertex,
        reshape(tb.range,2,nn),reshape(tb.child,2,nn),
        reshape(tb.lo,nr,nn),reshape(tb.hi,nr,nn),reshape(tb.b,nr,nn),tb.e)
end

function cem_build_node!(tb::CEMTreeBuilder,first,last)
    x,g,perm = tb.x,tb.g,tb.perm
    N = size(x,1)
    nr = max(N - 1,1)
    id = length(tb.e) + 1
    append!(tb.range,(first,last))
    append!(tb.child,(0,0))
    cols = @view perm[first:last]
    n = length(cols)
    #bounding box and least squares affine fit g ≈ a + b⋅x (reduced coordinates)
    lo,hi = fill(Inf,nr),fill(-Inf,nr)
    x̄ = zeros(nr)
    ḡ = 0.0
    for c in cols
        for i in 1:nr
            xi = N == 1 ? 0.0 : x[i,c]
            lo[i] = min(lo[i],xi)
            hi[i] = max(hi[i],xi)
            x̄[i] += xi
        end
        ḡ += g[c]
    end
    x̄ ./= n
    ḡ /= n
    S = zeros(nr,nr)
    t = zeros(nr)
    if N > 1
        for c in cols
            for i in 1:nr
                dxi = x[i,c] - x̄[i]
                t[i] += dxi*(g[c] - ḡ)
                for j in 1:nr
                    S[i,j] += dxi*(x[j,c] - x̄[j])
                end
            end
        end
    end
    trS = tr(S)
    for i in 1:nr
        S[i,i] += 1e-10*trS + 1e-300 #ridge: nodes can be degenerate in some directions
    end
    b = trS > 0 ? S \ t : zeros(nr)
    all(isfinite,b) || (b .= 0)
    e = Inf
    for c in cols
        ec = g[c]
        for i in 1:(N - 1)
            ec -= b[i]*x[i,c]
        end
        e = min(e,ec)
    end
    append!(tb.lo,lo)
    append!(tb.hi,hi)
    append!(tb.b,b)
    push!(tb.e,e)
    n <= tb.leafsize && return id
    #split: first by aggregate state, then along the widest dimension
    s1 = tb.state[cols[1]]
    if any(c -> tb.state[c] != s1,cols)
        sort!(cols,by = c -> (tb.state[c] != s1))
        mid = first + count(c -> tb.state[c] == s1,cols) - 1
    else
        w,dim = findmax(hi .- lo)
        iszero(w) && return id
        sort!(cols,by = c -> x[dim,c])
        mid = (first + last) ÷ 2
    end
    left = cem_build_node!(tb,first,mid)
    right = cem_build_node!(tb,mid + 1,last)
    tb.child[2id - 1] = left
    tb.child[2id] = right
    return id
end

#rigorous lower bound of gₖ - xₖ⋅μ on a node. μr = μ[1:N-1] .- μ[N]
@inline function cem_node_bound(env::CEMEnvelope,nd,μr,μN)
    @inbounds begin
        lb = env.node_e[nd] - μN
        for i in eachindex(μr)
            v = env.node_b[i,nd] - μr[i]
            lb += min(env.node_lo[i,nd]*v,env.node_hi[i,nd]*v)
        end
    end
    return lb
end

@inline function cem_reduced_cost(env::CEMEnvelope,c,μ)
    @inbounds begin
        r = env.g[c]
        for i in eachindex(μ)
            r -= env.x[i,c]*μ[i]
        end
    end
    return r
end

## LP workspace

struct CEMWorkspace
    B::Vector{Int32}
    Binv::Matrix{Float64}
    Bm::Matrix{Float64}
    λ::Vector{Float64}
    μ::Vector{Float64}
    μr::Vector{Float64}
    d::Vector{Float64}
    gB::Vector{Float64}
    stack::Vector{Int32}
    bstack::Vector{Float64}
    cand::Vector{Int32}
    kbuf::Vector{Int32}
end

function CEMWorkspace(env::CEMEnvelope)
    N = length(env)
    return CEMWorkspace(zeros(Int32,N),zeros(N,N),zeros(N,N),zeros(N),zeros(N),zeros(max(N - 1,0)),zeros(N),zeros(N),Int32[],Float64[],Int32[],zeros(Int32,N))
end

#most negative reduced cost over all columns (or the first one in tree order if `bland`), 0 if the basis is optimal
function cem_price!(ws::CEMWorkspace,env::CEMEnvelope,bland::Bool,::Nothing)
    μ,μr,stack,bstack = ws.μ,ws.μr,ws.stack,ws.bstack
    N = length(μ)
    μN = μ[N]
    for i in 1:(N - 1)
        μr[i] = μ[i] - μN
    end
    best_r = -CEM_RTOL
    best_c = 0
    empty!(stack)
    empty!(bstack)
    push!(stack,1)
    push!(bstack,cem_node_bound(env,1,μr,μN))
    @inbounds while !isempty(stack)
        nd = pop!(stack)
        pop!(bstack) >= best_r && continue #best_r may have improved since the node was pushed
        left,right = env.node_child[1,nd],env.node_child[2,nd]
        if left == 0
            for c in env.node_range[1,nd]:env.node_range[2,nd]
                r = cem_reduced_cost(env,c,μ)
                if r < best_r
                    best_r,best_c = r,c
                    bland && return best_c
                end
            end
        else
            bl,br = cem_node_bound(env,left,μr,μN),cem_node_bound(env,right,μr,μN)
            if bland || bl < br #visit the leftmost (bland) or most promising child first
                br < best_r && (push!(stack,right); push!(bstack,br))
                bl < best_r && (push!(stack,left); push!(bstack,bl))
            else
                bl < best_r && (push!(stack,left); push!(bstack,bl))
                br < best_r && (push!(stack,right); push!(bstack,br))
            end
        end
    end
    return best_c
end

#pricing restricted to the candidate columns
function cem_price!(ws::CEMWorkspace,env::CEMEnvelope,bland::Bool,cand::Vector{Int32})
    best_r = -CEM_RTOL
    best_c = 0
    for c in cand
        r = cem_reduced_cost(env,c,ws.μ)
        if r < best_r
            best_r,best_c = r,c
            bland && return best_c
        end
    end
    return best_c
end

function cem_refactor!(ws::CEMWorkspace,env::CEMEnvelope,z)
    for j in eachindex(ws.B)
        ws.Bm[:,j] .= @view(env.x[:,ws.B[j]])
    end
    F = lu!(ws.Bm,check = false)
    if !issuccess(F) || ws.B == env.vertex
        #restart from the pure component basis (B = I, always feasible)
        ws.B .= env.vertex
        ws.Binv .= 0
        for i in eachindex(ws.B)
            ws.Binv[i,i] = 1
        end
    else
        ws.Binv .= inv(F)
    end
    mul!(ws.λ,ws.Binv,z)
    ws.λ .= max.(ws.λ,0.0)
    return ws
end

#basis of pure component vertices: always feasible, λ = z
function cem_reset!(ws::CEMWorkspace,env::CEMEnvelope,z)
    ws.B .= env.vertex
    return cem_refactor!(ws,env,z)
end

const CEM_PIVTOL = 1e-9

#replace column r of the basis by c, given d = B⁻¹xc
function cem_pivot!(ws::CEMWorkspace,r,c)
    Binv,d = ws.Binv,ws.d
    N = length(d)
    @inbounds begin
        piv = d[r]
        for j in 1:N
            Binv[r,j] /= piv
        end
        for j in 1:N
            Brj = Binv[r,j]
            for i in 1:N
                i == r && continue
                Binv[i,j] -= d[i]*Brj
            end
        end
    end
    ws.B[r] = c
    return ws
end

#warm start for a feed equal to the composition of column c: put c in the basis (λ = eᵣ, feasible).
function cem_crash!(ws::CEMWorkspace,env::CEMEnvelope,c,z)
    ws.λ .= 0
    r = findfirst(isequal(c),ws.B)
    if r !== nothing
        ws.λ[r] = 1
        return ws
    end
    cem_mul!(ws.d,ws.Binv,env.x,c)
    _,r = findmax(abs,ws.d)
    if abs(ws.d[r]) < 1e-8
        return cem_reset!(ws,env,z)
    end
    cem_pivot!(ws,r,c)
    ws.λ[r] = 1
    return ws
end

#d = Binv*x[:,c] (small matrices: avoids BLAS overhead)
@inline function cem_mul!(d,Binv,x,c)
    N = length(d)
    @inbounds for i in 1:N
        s = 0.0
        for j in 1:N
            s += Binv[i,j]*x[j,c]
        end
        d[i] = s
    end
    return d
end

#revised primal simplex, starting from the (feasible) basis in ws.
#cand = nothing: all columns (tree pricing). Otherwise only the candidate columns (warm start, stops silently).
function cem_lp!(ws::CEMWorkspace,env::CEMEnvelope,z,cand = nothing)
    B,Binv,λ,μ,d,gB = ws.B,ws.Binv,ws.λ,ws.μ,ws.d,ws.gB
    N = length(B)
    degenerate = 0
    maxiter = cand === nothing ? 100_000 : 100N
    for iter in 1:maxiter
        iter % 64 == 0 && cem_refactor!(ws,env,z)
        @inbounds for i in 1:N
            s = 0.0
            for j in 1:N
                s += Binv[j,i]*env.g[B[j]]
            end
            μ[i] = s
        end
        c = cem_price!(ws,env,degenerate > 50,cand)
        c == 0 && return ws
        cem_mul!(d,Binv,env.x,c)
        r = 0
        θ = Inf
        for i in 1:N
            if d[i] > CEM_PIVTOL
                t = λ[i]/d[i]
                if t < θ || (t == θ && B[i] < B[r])
                    θ,r = t,i
                end
            end
        end
        r == 0 && throw(error("CEM: unbounded linear program (this should not happen)"))
        degenerate = θ <= 1e-14 ? degenerate + 1 : 0
        for i in 1:N
            λ[i] = max(λ[i] - θ*d[i],0.0)
        end
        λ[r] = θ
        cem_pivot!(ws,r,c)
    end
    cand === nothing && throw(error("CEM: maximum number of simplex iterations reached"))
    return ws
end

#index of a lattice point (in the order of cem_lattice)
function cem_lattice_rank(k,δ)
    N = length(k)
    R = δ
    idx = 1
    for i in 1:(N - 1)
        d = R - k[i]
        d > 0 && (idx += binomial(d - 1 + N - i,N - i))
        R -= k[i]
    end
    return idx
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
            cem_add_point_columns!(cand,env,cem_lattice_rank(kbuf,env.δ))
        end
    end
    return cand
end

#candidates for the local LP of lattice point k: neighbourhood of k and of the current basis
function cem_local_candidates!(ws::CEMWorkspace,env::CEMEnvelope,k)
    cand = ws.cand
    empty!(cand)
    cem_add_neighbourhood!(cand,env,k,ws.kbuf)
    for c in ws.B
        cem_add_neighbourhood!(cand,env,env.point[c],ws.kbuf)
    end
    return cand
end

#tangent plane (μᵢ/RT, in the reference of the model) of the current basis
function cem_tangent_plane(ws::CEMWorkspace,env::CEMEnvelope)
    for j in eachindex(ws.B)
        ws.gB[j] = env.g[ws.B[j]]
    end
    return transpose(ws.Binv)*ws.gB .+ env.gpure
end

## phases

#columns a and b can belong to the same phase (aggregate state)
function cem_compatible(env::CEMEnvelope,a,b)
    sa,sb = env.state[a],env.state[b]
    (sa == CEM_SOLID || sb == CEM_SOLID) && return a == b
    sa == sb && return true
    return isapprox(env.volume[a],env.volume[b],rtol = CEM_VOL_RTOL) #near critical: the roots merge
end

#stability (single phase) of lattice points, from a diagram or lazily from the LP of each point
struct CEMStability
    env::CEMEnvelope
    point_facet::Union{Nothing,Vector{Int32}}
    ws::CEMWorkspace
    memo::Dict{Int,Bool}
end

CEMStability(env::CEMEnvelope,point_facet = nothing) = CEMStability(env,point_facet,CEMWorkspace(env),Dict{Int,Bool}())

function (stable::CEMStability)(k)
    stable.point_facet !== nothing && return iszero(stable.point_facet[k])
    return get!(stable.memo,k) do
        env,ws = stable.env,stable.ws
        z = env.x[:,env.point_column[k]]
        cem_reset!(ws,env,z)
        cem_lp!(ws,env,z)
        cem_numphases(env,ws) == 1
    end
end

#lattice point next to x (in lattice units, sum(x) = δ)
function cem_round_lattice!(k,x,δ)
    k .= floor.(Int32,x)
    for _ in 1:(δ - sum(k))
        _,i = findmax(i -> x[i] - k[i],eachindex(x))
        k[i] += 1
    end
    return k
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
        cem_round_lattice!(kq,x,env.δ)
        ok = cem_stable_compatible(env,cem_lattice_rank(kq,env.δ),a,stable)
        for jj in 1:N, ii in 1:N
            (ok || ii == jj || kq[jj] == 0) && continue
            kn .= kq
            kn[jj] -= 1
            kn[ii] += 1
            ok = cem_stable_compatible(env,cem_lattice_rank(kn,env.δ),a,stable)
        end
        ok || return false
    end
    return true
end

cem_stable_compatible(env,q,a,stable) = stable(q) && cem_compatible(env,env.point_column[q],a)

#cluster index of each column in cols
function cem_clusters(env::CEMEnvelope,cols,stable = nothing)
    n = length(cols)
    parent = collect(1:n)
    root(i) = (while parent[i] != i; i = parent[i]; end; i)
    for a in 1:n, b in (a + 1):n
        root(a) == root(b) && continue
        if cem_samephase(env,cols[a],cols[b],stable)
            parent[root(b)] = root(a)
        end
    end
    roots = unique(root.(1:n))
    return [findfirst(isequal(root(i)),roots) for i in 1:n],length(roots)
end

function cem_numphases(env::CEMEnvelope,ws::CEMWorkspace,stable = nothing)
    cols = ws.B[ws.λ .> CEM_λTOL]
    length(cols) <= 1 && return 1
    return last(cem_clusters(env,cols,stable))
end

#phases (sorted by volume, solids last) of the current basis with weights λ
function cem_phases(env::CEMEnvelope,fluid,p,T,B,λ,stable = nothing)
    support = findall(>(CEM_λTOL),λ)
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
        comps[l] .+= w[j] .* @view(env.x[:,c])
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

function cem_flash_result(env::CEMEnvelope,fluid,p,T,z,ws::CEMWorkspace,stable = nothing)
    comps,β,vols,states = cem_phases(env,fluid,p,T,ws.B,ws.λ,stable)
    g = dot(ws.λ,@view(env.g[ws.B])) + dot(z,env.gpure)
    vapour_idx = -1
    for i in eachindex(states)
        states[i] == CEM_VAPOUR && (vapour_idx = i) #largest volume
    end
    return FlashResult(comps,β,vols,FlashData(p,T,g,vapour_idx))
end

## phase diagram

"""
    CEMPhaseDiagram(;δ = nothing, equilibrium = :auto, polish = nothing, threaded = true)

Method to calculate the complete phase diagram at constant pressure and temperature with the convex envelope method (CEM) [1].

The reduced Gibbs energy of every aggregate state (liquid, vapour and, for `CompositeModel` with a `SolidHfus` solid model, pure solids) is evaluated on a lattice of the composition simplex with step `1/δ`.
The lower convex envelope of those points is then determined for every lattice point. Its facets that connect points of different phases are the multiphase regions (tie-lines, tie-triangles, …).
Instead of building the complete convex hull, the envelope is obtained by solving a linear program for every lattice point, warm-started from the previous point. The tangent plane distances needed by the simplex method are evaluated with a tree that discards most of the grid with a rigorous bound.

## Keyword arguments
- `δ`: number of divisions of the composition axes. By default, `256`, `128` and `64` for binary, ternary and quaternary mixtures (as in [1]), and `32` for more components. The number of grid points is `binomial(δ + N - 1, N - 1)`.
- `equilibrium`: `:auto` considers liquid and vapour phases, `:lle` only liquid phases.
- `polish`: `nothing` (default) or a `TPFlashMethod` (e.g. `MichelsenTPFlash()` or `MultiPhaseTPFlash()`). If given, the phases of every multiphase facet are refined with that method, starting from the CEM phases. The refined result is only accepted if it lowers the Gibbs energy. Facets with solid phases are not refined.
- `threaded`: evaluate the grid and the envelope in parallel (uses `Threads.nthreads()` tasks).

## References
1. Göttl, Q., Rosen, N., & Burger, J. (2025). Convex envelope method for T, p flash calculations for mixtures with an arbitrary number of components and arbitrary aggregate states. Computers & Chemical Engineering, 203, 109326. [doi:10.1016/j.compchemeng.2025.109326](https://doi.org/10.1016/j.compchemeng.2025.109326)
"""
Base.@kwdef struct CEMPhaseDiagram <: ThermodynamicMethod
    δ::Union{Int,Nothing} = nothing
    equilibrium::Symbol = :auto
    polish::Union{Nothing,TPFlashMethod} = nothing
    threaded::Bool = true
end

"""
    CEMFacet

Multiphase facet of the convex envelope. Fields:
- `columns`: grid columns (indices into the envelope) of the facet vertices.
- `μ`: tangent plane of the facet, `μᵢ/RT` in the reference state of the model.
- `compositions`, `volumes`, `states`: phases in equilibrium (sorted by volume, solids last, with `NaN` volume). Without polishing, each phase composition is the mean of its facet vertices.
"""
struct CEMFacet
    columns::Vector{Int32}
    μ::Vector{Float64}
    compositions::Vector{Vector{Float64}}
    volumes::Vector{Float64}
    states::Vector{Symbol}
end

"""
    CEMDiagram

Result of [`phase_diagram`](@ref). Fields:
- `components`, `p`, `T`: system and conditions.
- `envelope`: grid and search tree. `envelope.lattice ./ envelope.δ` are the grid compositions (one column per point).
- `facets`: vector of [`Clapeyron.CEMFacet`](@ref), the multiphase facets of the convex envelope.
- `point_facet`: for each grid point, the index of the multiphase facet that contains it (`0` if the point is a stable single phase).
- `equilibrium`: equilibrium type used to build the diagram.

A `CEMDiagram` can be used for fast flash calculations at the same conditions via `CEMFlash(diagram = diagram)`.
"""
struct CEMDiagram
    components::Vector{String}
    p::Float64
    T::Float64
    equilibrium::Symbol
    envelope::CEMEnvelope
    facets::Vector{CEMFacet}
    point_facet::Vector{Int32}
end

Base.show(io::IO,env::CEMEnvelope) = print(io,"CEMEnvelope(δ = ",env.δ,", ",size(env.lattice,2)," grid points, ",length(env.g)," columns)")
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
Returns a [`Clapeyron.CEMDiagram`](@ref), containing all multiphase regions (as facets of the convex envelope of the Gibbs energy) and the phase classification of every grid point.

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
    δ = something(method.δ,cem_default_δ(N))
    env = cem_envelope(model,p,T,δ,method.equilibrium,method.threaded)
    fluid = cem_models(model,p,T,fill(1/N,N),method.equilibrium)[2]
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
    ws = CEMWorkspace(env)
    w = fill(1/N,N)
    facets = map(facet_cols) do cols
        ws.B .= cols
        z = env.x[:,cols]*w
        cem_refactor!(ws,env,z)
        μ = cem_tangent_plane(ws,env)
        comps,_,vols,states = cem_phases(env,fluid,p,T,cols,w,stable)
        f = CEMFacet(cols,μ,comps,vols,cem_state_symbol.(states))
        method.polish === nothing ? f : cem_polish_facet(model,p,T,z,f,env,fluid,ws,stable,method.polish)
    end
    return CEMDiagram(component_list(model),p,T,method.equilibrium,env,facets,point_facet)
end

#second pass over the multiphase facets found by the walk: vertices that are the same phase but not lattice neighbours
#(see cem_samephase) are merged, points and facets that become single phase are removed.
function cem_merge_phases(env::CEMEnvelope,facet_cols,point_facet)
    N = length(env)
    stable = CEMStability(env,copy(point_facet))
    labels = [first(cem_clusters(env,cols,stable)) for cols in facet_cols]
    Binv = [inv(env.x[:,cols]) for cols in facet_cols]
    λ = zeros(N)
    newpf = zeros(Int32,length(point_facet))
    used = falses(length(facet_cols))
    for k in eachindex(point_facet)
        f = point_facet[k]
        f == 0 && continue
        mul!(λ,Binv[f],@view(env.x[:,env.point_column[k]]))
        l = labels[f]
        if length(unique(l[i] for i in 1:N if λ[i] > CEM_λTOL)) > 1
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
    ws = CEMWorkspace(env)
    keys = Dict{Vector{Int32},Int}()
    list = Vector{Int32}[]
    proven = Set{Vector{Int32}}()
    resolved = falses(size(env.lattice,2))
    key = similar(ws.B)
    z = zeros(length(env))
    λ = ws.λ
    current_proven = false
    for (j,k) in enumerate(idx)
        resolved[k] && continue
        c = env.point_column[k]
        z .= @view(env.x[:,c])
        inside = false
        if current_proven
            mul!(λ,ws.Binv,z)
            inside = all(>=(-1e-12),λ)
            inside && (λ .= max.(λ,0.0))
        end
        if !inside
            if j == 1 || !current_proven
                cem_reset!(ws,env,z)
            else
                cand = cem_local_candidates!(ws,env,k) #before the crash: neighbourhood of the previous basis
                cem_crash!(ws,env,c,z)
                cem_lp!(ws,env,z,cand)
            end
            sort!(key .= ws.B)
            if !(key in proven)
                cem_lp!(ws,env,z)
                sort!(key .= ws.B)
                push!(proven,copy(key))
            end
            for cB in ws.B
                resolved[env.point[cB]] = true
            end
            current_proven = true
        end
        if cem_numphases(env,ws) > 1
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

function cem_polish_facet(model,p,T,z,f::CEMFacet,env,fluid,ws,stable,method)
    ws.λ .= 1/length(z)
    res = cem_flash_result(env,fluid,p,T,z,ws,stable)
    res2 = cem_polish(model,p,T,z,res,method)
    res2 === res && return f
    np = numphases(res2,true)
    np == length(f.compositions) || return f
    idx = [i for i in 1:numphases(res2) if is_active_phase(res2,i)]
    vapour_idx = res2.data.vapour_idx
    states = [i == vapour_idx ? :vapour : :liquid for i in idx]
    return CEMFacet(f.columns,f.μ,[collect(res2.compositions[i]) for i in idx],res2.volumes[idx],states)
end

include("tp_flash/CEMtp_flash.jl")

export phase_diagram, CEMPhaseDiagram
