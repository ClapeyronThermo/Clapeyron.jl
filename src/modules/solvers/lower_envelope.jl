#=
Lower convex envelope of a point cloud (xₖ,gₖ), xₖ in the unit simplex: for a given z,

    min Σₖ λₖ gₖ   s.t.  Σₖ λₖ xₖ = z,  λ ≥ 0

solved with a revised primal simplex method on an N×N basis. The unit vectors eᵢ must be points of the cloud
(the basis of the lowest of them is always feasible).
Pricing (the most negative reduced cost gₖ - xₖ⋅μ) is accelerated with a kd-tree that stores, for every node,
an affine under-estimator of g. This gives a rigorous lower bound of the reduced cost on the node.
=#

const ENVELOPE_RTOL = 1e-10 #reduced cost tolerance
const ENVELOPE_PIVTOL = 1e-9

#all points k (sum(k) = δ) of the simplex lattice, in lexicographically decreasing order
function simplex_lattice(N,δ)
    K = Matrix{Int32}(undef,N,binomial(δ + N - 1,N - 1))
    k = zeros(Int32,N)
    _simplex_lattice!(K,k,1,δ,Ref(0))
    return K
end

function _simplex_lattice!(K,k,i,rem,j)
    if i == length(k)
        k[i] = rem
        j[] += 1
        K[:,j[]] .= k
        return nothing
    end
    for ki in rem:-1:0
        k[i] = ki
        _simplex_lattice!(K,k,i + 1,rem - ki,j)
    end
    return nothing
end

#index of a lattice point (in the order of simplex_lattice)
function simplex_lattice_rank(k,δ)
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

#lattice point next to x (in lattice units, sum(x) = δ)
function simplex_lattice_round!(k,x,δ)
    k .= floor.(Int32,x)
    for _ in 1:(δ - sum(k))
        _,i = findmax(i -> x[i] - k[i],eachindex(x))
        k[i] += 1
    end
    return k
end

#points (columns) are stored in tree order
struct LowerEnvelope
    x::Matrix{Float64} #N × ncolumns
    g::Vector{Float64}
    vertex::Vector{Int32} #lowest column at each unit vector
    node_range::Matrix{Int32} #2 × nnodes
    node_child::Matrix{Int32} #2 × nnodes, 0 for leaves
    node_lo::Matrix{Float64} #(N-1) × nnodes, bounding box in reduced coordinates (x₁ … xₙ₋₁)
    node_hi::Matrix{Float64}
    node_b::Matrix{Float64} #slope of the affine under-estimator
    node_e::Vector{Float64} #minimum of g - b⋅x on the node
end

Base.length(env::LowerEnvelope) = size(env.x,1)

struct EnvelopeTreeBuilder
    x::Matrix{Float64}
    g::Vector{Float64}
    group::Vector{Int8}
    perm::Vector{Int32}
    range::Vector{Int32}
    child::Vector{Int32}
    lo::Vector{Float64}
    hi::Vector{Float64}
    b::Vector{Float64}
    e::Vector{Float64}
    leafsize::Int
end

#returns the envelope and the permutation of the columns of X (tree order). Points of different groups are split first.
function lower_envelope(X,G,group)
    N,nc = size(X)
    tb = EnvelopeTreeBuilder(X,G,group,collect(Int32,1:nc),Int32[],Int32[],Float64[],Float64[],Float64[],Float64[],16)
    _envelope_node!(tb,1,nc)
    perm = tb.perm
    x,g = X[:,perm],G[perm]
    vertex = [Int32(argmin(c -> x[i,c] == 1 ? g[c] : Inf,1:nc)) for i in 1:N]
    all(i -> x[i,vertex[i]] == 1,1:N) || throw(ArgumentError("lower_envelope: the unit vectors must be points of the cloud"))
    nn,nr = length(tb.e),N - 1
    env = LowerEnvelope(x,g,vertex,reshape(tb.range,2,nn),reshape(tb.child,2,nn),
        reshape(tb.lo,nr,nn),reshape(tb.hi,nr,nn),reshape(tb.b,nr,nn),tb.e)
    return env,perm
end

function _envelope_node!(tb::EnvelopeTreeBuilder,first,last)
    nr = size(tb.x,1) - 1
    id = length(tb.e) + 1
    append!(tb.range,(first,last))
    append!(tb.child,(0,0))
    cols = @view tb.perm[first:last]
    xr,gc = tb.x[1:nr,cols],tb.g[cols]
    lo,hi = vec(minimum(xr,dims = 2)),vec(maximum(xr,dims = 2))
    #least squares fit g ≈ a + b⋅x. any b gives a rigorous bound, e = min(g - b⋅x)
    b = (qr([ones(length(cols)) xr'],ColumnNorm()) \ gc)[2:end]
    all(isfinite,b) || (b .= 0)
    append!(tb.lo,lo)
    append!(tb.hi,hi)
    append!(tb.b,b)
    push!(tb.e,minimum(gc .- xr'*b))
    length(cols) <= tb.leafsize && return id
    #split: first by group, then along the widest dimension
    s1 = tb.group[cols[1]]
    if any(c -> tb.group[c] != s1,cols)
        sort!(cols,by = c -> (tb.group[c] != s1))
        mid = first + count(c -> tb.group[c] == s1,cols) - 1
    else
        w,dim = findmax(hi .- lo)
        iszero(w) && return id
        sort!(cols,by = c -> tb.x[dim,c])
        mid = (first + last) ÷ 2
    end
    left = _envelope_node!(tb,first,mid)
    right = _envelope_node!(tb,mid + 1,last)
    tb.child[2id - 1] = left
    tb.child[2id] = right
    return id
end

#rigorous lower bound of gₖ - xₖ⋅μ on a node. μr = μ[1:N-1] .- μ[N]
@inline function _envelope_node_bound(env::LowerEnvelope,nd,μr,μN)
    @inbounds begin
        lb = env.node_e[nd] - μN
        for i in eachindex(μr)
            v = env.node_b[i,nd] - μr[i]
            lb += min(env.node_lo[i,nd]*v,env.node_hi[i,nd]*v)
        end
    end
    return lb
end

@inline function _envelope_reduced_cost(env::LowerEnvelope,c,μ)
    @inbounds begin
        r = env.g[c]
        for i in eachindex(μ)
            r -= env.x[i,c]*μ[i]
        end
    end
    return r
end

struct EnvelopeWorkspace
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
end

function EnvelopeWorkspace(env::LowerEnvelope)
    N = length(env)
    return EnvelopeWorkspace(zeros(Int32,N),zeros(N,N),zeros(N,N),zeros(N),zeros(N),zeros(N - 1),zeros(N),zeros(N),Int32[],Float64[])
end

#most negative reduced cost over all columns (or the first one in tree order if `bland`), 0 if the basis is optimal
function _envelope_price!(ws::EnvelopeWorkspace,env::LowerEnvelope,bland::Bool,::Nothing)
    μ,μr,stack,bstack = ws.μ,ws.μr,ws.stack,ws.bstack
    N = length(μ)
    μN = μ[N]
    for i in 1:(N - 1)
        μr[i] = μ[i] - μN
    end
    best_r = -ENVELOPE_RTOL
    best_c = 0
    empty!(stack)
    empty!(bstack)
    push!(stack,1)
    push!(bstack,_envelope_node_bound(env,1,μr,μN))
    @inbounds while !isempty(stack)
        nd = pop!(stack)
        pop!(bstack) >= best_r && continue #best_r may have improved since the node was pushed
        left,right = env.node_child[1,nd],env.node_child[2,nd]
        if left == 0
            for c in env.node_range[1,nd]:env.node_range[2,nd]
                r = _envelope_reduced_cost(env,c,μ)
                if r < best_r
                    best_r,best_c = r,c
                    bland && return best_c
                end
            end
        else
            bl,br = _envelope_node_bound(env,left,μr,μN),_envelope_node_bound(env,right,μr,μN)
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
function _envelope_price!(ws::EnvelopeWorkspace,env::LowerEnvelope,bland::Bool,cand::Vector{Int32})
    best_r = -ENVELOPE_RTOL
    best_c = 0
    for c in cand
        r = _envelope_reduced_cost(env,c,ws.μ)
        if r < best_r
            best_r,best_c = r,c
            bland && return best_c
        end
    end
    return best_c
end

function envelope_refactor!(ws::EnvelopeWorkspace,env::LowerEnvelope,z)
    for j in eachindex(ws.B)
        ws.Bm[:,j] .= @view(env.x[:,ws.B[j]])
    end
    F = lu!(ws.Bm,check = false)
    if !issuccess(F) || ws.B == env.vertex
        #restart from the unit vector basis (B = I, always feasible)
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

#basis of unit vectors: always feasible, λ = z
function envelope_reset!(ws::EnvelopeWorkspace,env::LowerEnvelope,z)
    ws.B .= env.vertex
    return envelope_refactor!(ws,env,z)
end

#replace column r of the basis by c, given d = B⁻¹xc
function _envelope_pivot!(ws::EnvelopeWorkspace,r,c)
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

#warm start for a feed equal to column c: put c in the basis (λ = eᵣ, feasible).
function envelope_crash!(ws::EnvelopeWorkspace,env::LowerEnvelope,c,z)
    ws.λ .= 0
    r = findfirst(isequal(c),ws.B)
    if r !== nothing
        ws.λ[r] = 1
        return ws
    end
    _envelope_mul!(ws.d,ws.Binv,env.x,c)
    _,r = findmax(abs,ws.d)
    if abs(ws.d[r]) < 1e-8
        return envelope_reset!(ws,env,z)
    end
    _envelope_pivot!(ws,r,c)
    ws.λ[r] = 1
    return ws
end

#d = Binv*x[:,c] (small matrices: avoids BLAS overhead)
@inline function _envelope_mul!(d,Binv,x,c)
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
function envelope_lp!(ws::EnvelopeWorkspace,env::LowerEnvelope,z,cand = nothing)
    B,Binv,λ,μ,d = ws.B,ws.Binv,ws.λ,ws.μ,ws.d
    N = length(B)
    degenerate = 0
    maxiter = cand === nothing ? 100_000 : 100N
    for iter in 1:maxiter
        iter % 64 == 0 && envelope_refactor!(ws,env,z)
        @inbounds for i in 1:N
            s = 0.0
            for j in 1:N
                s += Binv[j,i]*env.g[B[j]]
            end
            μ[i] = s
        end
        c = _envelope_price!(ws,env,degenerate > 50,cand)
        c == 0 && return ws
        _envelope_mul!(d,Binv,env.x,c)
        r = 0
        θ = Inf
        for i in 1:N
            if d[i] > ENVELOPE_PIVTOL
                t = λ[i]/d[i]
                if t < θ || (t == θ && B[i] < B[r])
                    θ,r = t,i
                end
            end
        end
        r == 0 && error("lower envelope: unbounded linear program (this should not happen)")
        degenerate = θ <= 1e-14 ? degenerate + 1 : 0
        for i in 1:N
            λ[i] = max(λ[i] - θ*d[i],0.0)
        end
        λ[r] = θ
        _envelope_pivot!(ws,r,c)
    end
    cand === nothing && error("lower envelope: maximum number of simplex iterations reached")
    return ws
end

#tangent plane μ (Σ μᵢxᵢ = g on the facet) of the current basis
function envelope_tangent_plane(ws::EnvelopeWorkspace,env::LowerEnvelope)
    for j in eachindex(ws.B)
        ws.gB[j] = env.g[ws.B[j]]
    end
    return transpose(ws.Binv)*ws.gB
end
