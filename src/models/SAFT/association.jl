function a_assoc(model::EoSModel, V, T, z,data=nothing)
    _0 = zero(V+T+first(z))
    nn = assoc_pair_length(model)
    iszero(nn) && return _0
    isone(nn) && return a_assoc_exact_1(model,V,T,z,data)
    #_X,_Δ = @f(X_and_Δ,data)
    #return @f(a_assoc_impl,_X,_Δ)
    _X = @f(X,data)
    return @f(a_assoc_impl,_X)
end

"""
    assoc_shape(model::EoSModel)::Compressed4DMatrix{Int64,UnitRange{Int64}}

Returns a `Clapeyron.Compressed4DMatrix` that has the same shape as the association sites used by the model.
By default, it has the same shape as `model.params.bondvol`
## Example:

```julia-repl
julia> model = PCSAFT(["water"])
PCSAFT{BasicIdeal} with 1 component:
 "water"
Contains parameters: Mw, segment, sigma, epsilon, epsilon_assoc, bondvol

julia> Clapeyron.assoc_shape(model)
Clapeyron.Compressed4DMatrix{Int64, UnitRange{Int64}} with 1 entry:
 (1, 1) >=< (1, 2): 1 #component 1 at site 1 has association interaction with component 1 at site 2.
```
"""
assoc_shape(model::EoSModel) = assoc_shape(model.params.bondvol)
assoc_shape(param::AssocParam) = assoc_shape(param.values)
@inline function assoc_shape(mat::Compressed4DMatrix)
    l = length(mat.values)
    Compressed4DMatrix{Int64,UnitRange{Int64}}(1:l,mat.outer_indices,mat.inner_indices,mat.outer_size,mat.inner_size)
end

@inline function assoc_similar(model::EoSModel,::Type{𝕋}) where 𝕋
    assoc_similar(assoc_shape(model),𝕋)
end

assoc_similar(model::EoSModel) = assoc_similar(model,eltype(model))

"""
    assoc_pair_length(model::EoSModel)

Indicates the number of pair combinations between the different sites in an association model.
By default uses `length(assoc_shape(model).values)`

## Example:

```julia-repl
julia> model = PCSAFT(["water"])
PCSAFT{BasicIdeal} with 1 component:
 "water"
Contains parameters: Mw, segment, sigma, epsilon, epsilon_assoc, bondvol

julia> model.params.bondvol
AssocParam{Float64}["water"]) with 1 value:
("water", "e") >=< ("water", "H"): 0.034868

julia> Clapeyron.assoc_pair_length(model)
1
```
"""
@inline function assoc_pair_length(model::EoSModel)
    val = assoc_shape(model)
    return length(val.values)
end

"""
    assoc_strength(model::EoSModel,V,T,z,i,j,a,b,data = Clapeyron.data(Model,V,T,z))
    Δ(model::EoSModel,V,T,z,i,j,a,b,data = Clapeyron.data(Model,V,T,z))
Calculates the asssociation strength between component `i` at site `a` and component `j` at site `b`.

Any precomputed values can be passed along by calling `Clapeyron.data`.

## Example
```julia-repl
julia> model = PCSAFT(["water"])
PCSAFT{BasicIdeal} with 1 component:
 "water"
Contains parameters: Mw, segment, sigma, epsilon, epsilon_assoc, bondvol

julia> model.params.bondvol.values
Clapeyron.Compressed4DMatrix{Float64, Vector{Float64}} with 1 entry:
 (1, 1) >=< (1, 2): 0.034868

julia> Clapeyron.assoc_strength(model,2.5e-5,298.15,[1.0],1,1,1,2) #you can also use Clapeyron.Δ
1.293144062056963e-26

#PCSAFT precomputed data: (d,ζ₀,ζ₁,ζ₂,ζ₃,m̄)
julia> _data = Clapeyron.data(model,2.5e-5,298.15,[1.0])
([2.991688553098391e-10], 1.3440137996322956e28, 4.020870699566213e18, 1.2029192845380957e9, 0.3598759853853927, 1.0656)

julia> Clapeyron.Δ(model,2.5e-5,298.15,[1.0],1,1,1,2,_data)
1.293144062056963e-26
```
"""
function Δ end
const assoc_strength = Δ

"""
    Δ(model::EoSModel, V, T, z)
    Δ(model::EoSModel, V, T, z,data)
    assoc_strength(model::EoSModel, V, T, z)
    assoc_strength(model::EoSModel, V, T, z,data)

Returns a list of all combinations of non-zero association strength, calculated at V,T,z conditions. returns a `Clapeyron.Compressed4DMatrix`.
By default, it calls `assoc_similar(model,𝕋)` (where 𝕋 is the promoted type of all the arguments) and fills the list using `Δ(model,V,T,z,i,j,a,b,data)`

## Example
```julia-repl
julia> model = PCSAFT(["water"])
PCSAFT{BasicIdeal} with 1 component:
 "water"
Contains parameters: Mw, segment, sigma, epsilon, epsilon_assoc, bondvol

julia> model.params.bondvol
AssocParam{Float64}["water"]) with 1 value:
("water", "e") >=< ("water", "H"): 0.034868

julia> Clapeyron.assoc_strength(model,2.5e-5,298.15,[1.0],1,1,1,2) #you can also use Clapeyron.Δ
1.293144062056963e-26

#PCSAFT precomputed data: (d,ζ₀,ζ₁,ζ₂,ζ₃,m̄)
julia> _data = Clapeyron.data(model,2.5e-5,298.15,[1.0])
([2.991688553098391e-10], 1.3440137996322956e28, 4.020870699566213e18, 1.2029192845380957e9, 0.3598759853853927, 1.0656)

julia> Clapeyron.Δ(model,2.5e-5,298.15,[1.0],1,1,1,2,_data)
1.293144062056963e-26
```
"""
function Δ(model::EoSModel, V, T, z)
    Δout = assoc_similar(model,typeof(V+T+first(z)))
    for (idx,(i,j),(a,b)) in indices(Δout)
        Δout[idx] =@f(Δ,i,j,a,b)
    end
    return Δout
end

"""
    delta_assoc(model,V,T,z,data)

equivalent to `Δ(model,V,T,z)`, but additionally, performs runtime mixing of `Δ[i,j][a,b]` values if specified via `assoc_options`.

"""
function delta_assoc(model,V,T,z,data::M) where M
    if data === nothing
        delta = Δ(model,V,T,z)
    else
        delta = Δ(model,V,T,z,data)
    end
    options = assoc_options(model)
    combining = options.combining
    if combining in (:elliott_runtime,:esd_runtime)
        elliott_runtime_mix!(delta)
    end
    return delta
end

"""
    assoc_options(model::EoSModel)

Returns association options used in the association solver.

"""
@inline function assoc_options(model::EoSModel)
    return model.assoc_options
end

function Δ(model::EoSModel, V, T, z, data)
    Δout = assoc_similar(model,@f(Base.promote_eltype))
    Δout.values .= false
    for (idx,(i,j),(a,b)) in indices(Δout)
        Δout[idx] =@f(Δ,i,j,a,b,data)
    end
    return Δout
end

function issite(i::Int,a::Int,ij::Tuple{Int,Int},ab::Tuple{Int,Int})::Bool
    ia = (i,a)
    i1,i2 = ij
    a1,a2 = ab
    ia1 = (i1,a1)
    ia2 = (i2,a2)
    return (ia == ia1) | (ia == ia2)
end

function complement_index(i,ij)::Int
    i1,i2 = ij
    ifelse(i1 == i,i2,i1)::Int
end

function compute_index(idxs,i,a)::Int
    res::Int = idxs[i] + a - 1
    return res
end

function inverse_index(idxs,o)
    i = findfirst(>=(o-1),idxs)::Int
    a = o + 1 - idxs[i]
    return i,a
end

function assoc_site_matrix(model,V,T,z,data = nothing,delta = @f(delta_assoc,data))
    options = assoc_options(model)
    return dense_assoc_site_matrix(model,V,T,z,data,delta)
end

#this fills the zeros of the Δ vector with the corresponding mixing values
function elliott_runtime_mix!(Δ)
    _Δ = Δ.values
    for (idx1,(i1,i2),(a1,a2)) in indices(Δ)
        if i1 == i2
            i = i1
            Δi = _Δ[idx1]
            for (idx2,(j1,j2),(b1,b2)) in indices(Δ)
                if j1 == j2
                    j = j1
                    Δj = _Δ[idx2]
                    Δijab = sqrt(Δi*Δj)
                    if !iszero(Δijab)
                        Δij = Δ[i,j]
                        v_idx1 = validindex(Δij,a1,b2)
                        v_idx2 = validindex(Δij,a2,b1)
                        v_idx1 != 0 && iszero(_Δ[v_idx1]) && (_Δ[v_idx1] = Δijab)
                        v_idx2 != 0 && iszero(_Δ[v_idx2]) && (_Δ[v_idx2] = Δijab)
                    end
                end
            end
        end
    end
    return Δ
end

function dense_assoc_site_matrix(model,V,T,z,data=nothing,delta = @f(delta_assoc,data))
    sitesparam = getsites(model)
    _sites = sitesparam.n_sites
    p = _sites.p
    ρ = N_A/V
    _ii::Vector{Tuple{Int,Int}} = delta.outer_indices
    _aa::Vector{Tuple{Int,Int}} = delta.inner_indices
    _idx = 1:length(_ii)
    Δ = delta.values
    TT = eltype(Δ)
    _n = sitesparam.n_sites.v
    nn = length(_n)
    K  = zeros(TT,nn,nn)
    options = assoc_options(model)
    combining = options.combining
    runtime_combining = combining ∈ (:elliott_runtime,:esd_runtime)

    @inbounds for i ∈ 1:length(z) #for i ∈ comps
        sitesᵢ = 1:(p[i+1] - p[i]) #sites are normalized, with independent indices for each component
        for a ∈ sitesᵢ #for a ∈ sites(comps(i))
            ia = compute_index(p,i,a)
            for idx ∈ _idx #iterating for all sites
                ij = _ii[idx]
                ab = _aa[idx]
                if issite(i,a,ij,ab)
                    j = complement_index(i,ij)
                    b = complement_index(a,ab)
                    jb = compute_index(p,j,b)
                    njb = _n[jb]
                    zj = z[j]
                    if !iszero(zj)
                        K[ia,jb]  = ρ*njb*z[j]*Δ[idx]
                    end
                end
            end
        end
    end
    return K::Matrix{TT}
end

function X end
const assoc_fractions = X

"""
    assoc_fractions(model::EoSModel, V, T, z,data = nothing)

Returns the solution for the association site fractions. used internally by all models that require association.
The result is of type `PackedVectorsOfVectors.PackedVectorOfVectors`, with `length = length(model)`, and `x[i][a]` representing the empty fraction of the site `a` at component `i`
## Example:

```
julia> model = PCSAFT(["water","methanol","ethane"],assoc_options = AssocOptions(combining = :esd))
PCSAFT{BasicIdeal} with 3 components:
 "water"
 "methanol"
 "ethane"
Contains parameters: Mw, segment, sigma, epsilon, epsilon_assoc, bondvol

julia> x = Clapeyron.assoc_fractions(model,2.6e-5,300.15,[0.3,0.3,0.4]) #you can also use `Clapeyron.X`
3-element pack(::Vector{Vector{Float64}}):
 [0.041396427041509046, 0.041396427041509046]
 [0.018874664357682362, 0.018874664357682362]
 0-element view(::Vector{Float64}, 5:4) with eltype Float64
```
"""
function X(model::EoSModel, V, T, z,data = nothing)
    #we return X with derivative information
    nn = assoc_pair_length(model)
    isone(nn) && return X_exact1(model,V,T,z,data)
    X,Δ = X_and_Δ(model,V,T,z,data)
    return X
    #for some reason, this fails on infinite dilution derivatives
    #=
    if eltype(X.v) === eltype(Δ.values)
        return X
    end
    X̄ = X.v
    #K matrix with derivative information
    K = assoc_site_matrix(model,V,T,z,data,Δ)
    X̃ = similar(K,length(X̄))

    #=
    strategy to obtain general derivatives of nonbonded fractions with automatic differenciation:

    using Implicit AD, we can update X with a "perfect newton upgrade", with the derivative information added in the last update.
    it is equivalent to the method of Tan (2004), in the sense that we still need to solve a linear system of equations containing X.
    but this only requires to solve one linear system, as the derivatives are carried by the number type, instead of separated.

    =#
    mul!(X̃,K,X̄)
    K .*= -1
    for k in 1:size(K,1)
        K[k,k] -= (1 + X̃[k])/X̄[k]
    end
    X̃ .+= -1 ./ X̄ .+ 1

    F = Solvers.unsafe_LU!(K)
    ldiv!(F,X̃)
    X̃ .+= X̄
    return PackedVofV(X.p,X̃) =#
end

function X_and_Δ(model::EoSModel, V, T, z,data = nothing)
    nn = assoc_pair_length(model)
    isone(nn) && return X_and_Δ_exact1(model,V,T,z,data)
    options = assoc_options(model)::AssocOptions
    _Δ = delta_assoc(model,V,T,z,data)
    #K = assoc_site_matrix(model,primalval(V),T,primalval(z),data,primalval(_Δ))
    K = assoc_site_matrix(model,V,T,z,data,_Δ)
    sitesparam = getsites(model)
    idxs = sitesparam.n_sites.p
    compress = __maybe_compress(K)
    if compress
        J,J_to_K = compress_assoc_matrix(K)
        Ysol = assoc_matrix_solve(J,options)
        Xsol = Ysol[J_to_K]
    else
        Xsol = assoc_matrix_solve(K,options)
    end
    return PackedVofV(idxs,Xsol),_Δ
end

function assoc_matrix_solve(K::AbstractMatrix{T},options::AssocOptions = AssocOptions()) where T
    atol = options.atol
    rtol = options.rtol
    max_iters = options.max_iters
    α = options.dampingfactor
    implicit_ad = options.implicit_ad
    return assoc_matrix_solve(K, α, atol ,rtol, max_iters, implicit_ad)
end

function check_antidiagonal2(x::AbstractMatrix)
    size(x) == (2,2) || return false
    x11,x22 = x[1,1],x[2,2]
    x21,x12 = x[2,1],x[1,2]
    return iszero(x11) & iszero(x22) & (x12 >= 0) & (x21 >= 0)
end

function check_antidiagonal22(x::AbstractMatrix)
    size(x) == (4,4) || return false
    return check_antidiagonal2(@view(x[1:2,1:2])) & check_antidiagonal2(@view(x[1:2,3:4])) &
    check_antidiagonal2(@view(x[3:4,1:2])) & check_antidiagonal2(@view(x[3:4,3:4]))
end

function assoc_matrix_x0!(K,X)
    #(A*x .* x) + x - 1 = 0
    success = false
    init = false
    if size(K) == (1,1)
        #1-site association
        k = K[1,1]
        #axx + x - 1 = 0
        #-1 +- sqrt(1 + 4a)/2 = 0
        X[1] = 0.5*(-1 + sqrt(1 + 4k))
        success = true
        init = true
    elseif size(K) == (2,2) && all(!iszero,K)
        #this matrix type is normally generated via compression
        X_exact2_denseK!(K,X)
        init = true
        success = true
    elseif check_antidiagonal2(K)
        X_exact2!(K,X)
        init = true
        success = true
    elseif check_antidiagonal22(K)
    #nb-nb association with cross-association
    K11 = @view(K[1:2,1:2])
    K12 = @view(K[1:2,3:4])
    K21 = @view(K[3:4,1:2])
    K22 = @view(K[3:4,3:4])

    X_exact2!(K11,@view(X[1:2]))
    X_exact2!(K22,@view(X[3:4]))
    if (iszero(K12) & iszero(K21))
        #solve each association separately, if one of the diagonal association
        #submatrices is zero, then cross-association does not have any sense.
        success = true
    else 
        #general solution, takes longer to compile.
        #_,success = X_exact4!(K,X)
        #success || X_exact2!(K22,@view(X[3:4]))
        success = false
    end
    init = true
    else
        #TODO: add more exact expressions.
    end

    if !init
        Kmin,Kmax = nonzero_extrema(K) #look for 0 < Amin < Amax
        if Kmax > 1
            f = true/Kmin
        else
            f = true-Kmin
        end
        fill!(X,min(f,one(f)))
    end
    return X,success
end


#this function destroys KK and XX0
function __assoc_matrix_solve_static(::Val{N},KK::AbstractMatrix{T1},XX0::AbstractVector{T2}, α, atol ,rtol, max_iters) where {N,T1,T2}
    X0 = SVector{N,T2}(XX0)
    K = SMatrix{N,N,T1,N*N}(KK)
    Xsol = X0
    it_ss = (5*length(Xsol))
    converged = false
    for i in 1:it_ss
        kx = K*X0
        Xsol = α ./ (1 .+ kx) .+ (1 .- α) .* X0
        converged,finite = Solvers.convergence(Xsol,X0,atol,rtol)
        if converged
            if !finite
                Xsol = NaN .* Xsol
            end

            break
        end
        X0 = Xsol
    end

    if converged
        XX0 .= Xsol
        return XX0
    end

    H = KK
    g = XX0
    #TODO: for the next stable release, use MVector
    piv = zeros(Int,N)
    for i in (it_ss + 1):max_iters
        #@show Xsol
        KX = K*Xsol
        H .= 0
        H .= -K
        for k in 1:size(H,1)
            H[k,k] -= (1 + KX[k])/Xsol[k]
        end
        F = Solvers.unsafe_LU!(H,piv)
        g .= 1 ./ Xsol .- 1 .- KX #gradient
        ldiv!(F,g)
        ΔX = SVector{N,T2}(XX0)
        Xnewton = Xsol - ΔX
        Xss = 1 ./ (1 .+ KX)
        X0 = Xsol
        Xsol = ifelse.(0 .<= Xnewton .<= 1, Xnewton, Xss)
        converged,finite = Solvers.convergence(Xsol,X0,atol,rtol,false,Inf)
        #@show converged,finite
        if converged
            if !finite
                Xsol = NaN .* Xsol
            end
            XX0 .= Xsol
            break
        end
    end

    if !converged
        Xsol = NaN .* Xsol
    end
    XX0 .= Xsol
    return XX0
end

function assoc_matrix_solve(K::AbstractMatrix{T}, α, atol ,rtol, max_iters, implicit_ad) where T
    n = LinearAlgebra.checksquare(K) #size
    #initialization procedure:
    X0 = Vector{T}(undef,n)
    X0,success = assoc_matrix_x0!(K,X0)
    success && return X0
    #static versions to improve speed, but they cost too much compilation time
    #length(X0) == 3 && return __assoc_matrix_solve_static(Val{3}(), K, X0, α, atol ,rtol, max_iters)
    #length(X0) == 4 && return __assoc_matrix_solve_static(Val{4}(), K, X0, α, atol ,rtol, max_iters)
    #length(X0) == 5 && return __assoc_matrix_solve_static(Val{5}(), K, X0, α, atol ,rtol, max_iters)
    
    if implicit_ad && K[1] isa ForwardDiff.Dual
        K_primal = nested_pvalue(K) # solve on primalval
        
        Xsol = assoc_matrix_solve_general(K_primal, nested_pvalue.(X0), n, α, atol ,rtol, max_iters)
        return assoc_matrix_solve_ad(Xsol, K, K_primal) # implicit AD
    else
        #Propagate AD through solver
        return assoc_matrix_solve_general(K, X0, n, α, atol ,rtol, max_iters)
    end
end

"""
    J,idx = compress_assoc_matrix(K)

Given a matrix generated by `Clapeyron.dense_assoc_site_matrix`, `Clapeyron.compress_assoc_matrix` will check if there is an equivalent, smaller matrix J, that solves the same association system.
returns the modified matrix and a set of indices such as `assoc_matrix_solve(K) == assoc_matrix_solve(J)[idx]`

The compression strategy works mainly with association systems with donor-aceptor pairs containing the same amount of sites.

If compression fails, it returns unmodified indices and the original matrix.
"""
function compress_assoc_matrix(K)
    n = LinearAlgebra.checksquare(K)
    idx = zeros(Int,n)
    idx .= 1:n
    n_unique = n

    #step 1: check which rows are the same if we translate the row by one and pad with zeros.
    for i in 1:n
        !iszero(K[1,i]) && continue
        idx[i] != i && continue #already compressed
        Ki = @view K[2:end,i]
        for j in (i+1):n
            !iszero(K[end,j]) && continue
            Kj = @view K[1:end-1,j]
            if Ki == Kj
                n_unique -= 1
                idx[j] = i
            end
        end
    end
    #failed to compress
    if n_unique == n
        return K,idx
    end

    can_compress = true
    #step 2: check consistency, if we remove n rows, there it at least n-1 (translated) cols of pure zeros:
    for i in 1:n
        idx[i] != i && continue #already compressed
        for ii in 2:n
            Kiii = K[ii,i]
            if iszero(Kiii)
                for j in (i+1):n
                    if idx[j] == i
                       can_compress = iszero(K[ii-1,j])
                    end
                end
            end
        end
    end
    
    if !can_compress
        idx .= 1:n
        return K,idx
    end

    #step 3: perform compression
    J = similar(K,(n_unique,n_unique))
    ii = 0
    for i in 1:n
        if idx[i] == i
            ii += 1
            jj = 0
            for j in 1:n
                Kji = K[j,i]
                if !iszero(Kji)
                    jj += 1
                    J[jj,ii] = Kji
                end
            end
        end
    end

    #step 4: generate output indices
    i_unique = 0
    idx .*= -1
    for i in 1:n
        ix = idx[i]
        if ix < 0
            i_unique += 1
            idx[i] = i_unique
            for j in (i+1):n
                if idx[j] == ix
                    idx[j] = i_unique
                end
            end
        end
    end
    return J,idx
end

function __maybe_compress(K)
    n = LinearAlgebra.checksquare(K)
    #step 1: check which rows are the same if we translate the row by one and pad with zeros.
    for i in 1:n
        !iszero(K[1,i]) && continue
        Ki = @view K[2:end,i]
        for j in (i+1):n
            !iszero(K[end,j]) && continue
            Kj = @view K[1:end-1,j]
            if Ki == Kj
                return true
            end
        end
    end
   return false
end

function assoc_matrix_solve_general(K::AbstractMatrix{T}, X0, n, α, atol ,rtol, max_iters) where T
    Xsol = Vector{T}(undef,n)
    Xsol .= X0
    #=
    function to solve
    find vector x that satisfies:
    (A*x .* x) + x - 1 = 0
    solved by reformulating in succesive substitution:
    x .= 1 ./ (1 .+ A*x)

    #we perform a "partial multiplication". that is, we use the already calculated
    #values of the next Xi to calculate the current Xi. this seems to accelerate the convergence
    #by around 50% (check what the ass_matmul! function does)

    note that the damping is done inside the partial multiplication. if is done outside, it causes convergence problems.

    after a number of ss iterations are done, we use newton minimization.
    the code for the newton optimization is based on sgtpy: https://github.com/gustavochm/sgtpy/blob/336cb2a7581b22492914233e29062f5a364b47da/sgtpy/vrmie_pure/association_aux.py#L33-L57

    some notes:
    - the linear system is solved via LU decomposition, for that, we need to allocate one (1) Matrix{T} and one (1) Vector{Int}
    - gauss-seidel does not require an additional matrix allocation, but it is slow. (slower than SS)
    - julia 1.10 does not have a way to make LU non-allocating, but the code is simple, so it was added as the function unsafe_LU! in the Solvers module.
    =#
    fx(kx,x) =  α/(1+kx) + (1-α)*x
    function f_ss!(out,in)
        ass_matmul!(fx,out,K,in)
        return out
    end

    #successive substitution. 50 iters
    it_ss = (5*length(Xsol))
    converged = false
    for i in 1:it_ss
        f_ss!(Xsol,X0)
        converged,finite = Solvers.convergence(Xsol,X0,atol,rtol)
        if converged
            finite || (Xsol .= NaN)
            return Xsol
        end
        X0 .= Xsol
       # @show Xsol
    end
    if converged
        !finite && (Xsol .= NaN)
        return Xsol
    end
    H = Matrix{T}(undef,n,n)
    H .= 0
    piv = zeros(Int,n)
    dX = copy(Xsol)
    KX = copy(Xsol)
    for i in (it_ss + 1):max_iters
        #@show Xsol
        KX = mul!(KX,K,Xsol)
        H .= -K
        for k in 1:size(H,1)
            H[k,k] -= (1 + KX[k])/Xsol[k]
        end
        #F already contains H and the pivots, because we refreshed H, we need to refresh
        #the factorization too.
        F = Solvers.unsafe_LU!(H,piv)
        dX .= 1 ./ Xsol .- 1 .- KX #gradient
        ldiv!(F,dX) #we solve H/g, overwriting g
        X0 .= Xsol
        for k in 1:length(dX)
            Xk = Xsol[k]
            X0k = X0[k]
            dXk = dX[k]
            X0[k] = Xk
            X_newton = Xk - dXk
            if !(0 <= X_newton <= 1)
                Xsol[k] = 0.5*(Xk + X0k) #successive substitution step
            else
                Xsol[k] = X_newton #newton step
            end   
        end
        # Xsol .-= dX
        converged,finite = Solvers.convergence(Xsol,X0,atol,rtol,false,Inf)
        #@show converged,finite
        if converged
            finite || (Xsol .= NaN)
            return Xsol
        end
    end

    converged || (Xsol .= NaN)
    return Xsol
end

function assoc_matrix_solve_ad(Xsol::X, K::KT, K_primal::KP)::Vector{V2} where {V1,V2,X<:AbstractVector{V1},KT<:AbstractMatrix{V2},KP<:AbstractMatrix{V1}}
    N = Val(length(Xsol))
    f(X_::XX_,tups_::Tuple{KK_,Val{N_}}) where {V1_,V2_,N_,XX_<:AbstractVector{V1_},KK_<:AbstractMatrix{V2_}} = begin 
        TT = promote_type(V1_,V2_)
        _1 = one(TT)
        K_ = tups_[1]
        #itt = (@inbounds begin
        #            ki = @view K_[i,:]
        #            dot(ki,X_) * X_[i] + X_[i] - _1
        #        end
        #        for i in 1:N_) # Generator for inefficient code below
        #tmp = SVector{N_,TT}(itt)
        tmp = K_ * X_ # K * X, allocates initial Vector{TT} buffer=#
        tmp .*= X_ # (K * X) .* X
        tmp .+= X_ # (K * X) .* X .+ X
        tmp .-= _1 # (K * X) .* X .+ X .- 1
        tmp
    end# (1 + ∑_{jb} K⁽ⁱᵃʲᵇ⁾X⁽ⁱᵃʲᵇ⁾ )⁻¹ = Xⁱᵃ for all iₐ, but rearranged. 
    # f wrt X has polynomial form, which is easier (and more efficient) to differentiate compared to 1 ./ X
    return __gradients_for_root_finders(Xsol,(K,N),(K_primal,N),f) # implicit AD
end

#exact calculation of site non-bonded fraction when there is only one site

function X_exact1(model,V,T,z,data = nothing)
    xia,xjb,i,j,a,b,n,idxs,Δijab = _X_exact1(model,V,T,z,data)
    pack_X_exact1(xia,xjb,i,j,a,b,n,idxs)
end

function X_and_Δ_exact1(model,V,T,z,data = nothing)
    xia,xjb,i,j,a,b,n,idxs,Δijab = _X_exact1(model,V,T,z,data)
    XX = pack_X_exact1(primalval(xia),primalval(xjb),i,j,a,b,n,idxs)
    Δout = assoc_similar(model,@f(Base.promote_eltype))
    Δout.values[1] = Δijab
    return XX,Δout
end

function _X_exact1(model,V,T,z,data=nothing)
    κ = assoc_shape(model)
    i,j = κ.outer_indices[1]
    a,b = κ.inner_indices[1]
    if data === nothing
        _Δ = @f(Δ,i,j,a,b)
    else
        _Δ = @f(Δ,i,j,a,b,data)
    end
    _1 = one(eltype(_Δ))
    sitesparam = getsites(model)
    idxs = sitesparam.n_sites.p
    n = length(sitesparam.n_sites.v)
    ρ = N_A/V
    zi = z[i]
    zj = z[j]
    ni = sitesparam.n_sites[i]
    na = ni[a]
    nj = sitesparam.n_sites[j]
    nb = nj[b]
    ρ = N_A/V
    kia = na*zi*ρ*_Δ
    kjb = nb*zj*ρ*_Δ
    _a = kia
    _b = _1 - kia + kjb
    _c = -_1
    denom = _b + sqrt(_b*_b - 4*_a*_c)
    xia = -2*_c/denom
    xk_ia = kia*xia
    xjb = (1- xk_ia)/(1 - xk_ia*xk_ia)
    return xia,xjb,i,j,a,b,n,idxs,_Δ
end

function pack_X_exact1(xia,xjb,i,j,a,b,n,idxs)
    Xsol = fill(one(xia),n)
    _X = PackedVofV(idxs,Xsol)
    _X[j][b] = xjb
    _X[i][a] = xia
    return _X
end

#helper function to get the sites. in almost all cases, this is model.sites
#but SAFTgammaMie uses model.vrmodel.sites instead
"""
    getsites(model::EoSModel)

returns the `SiteParam` used in association calculations for the input `model`, if any. Fails if not available. to check if a model has sites, use `has_sites`

"""
getsites(model) = model.sites

function a_assoc_impl(model::EoSModel, V, T, z, X, Δ)
    #=
    todo: fix mixed derivatives at infinite dilution
    =#
    #=
    Implementation notes

    We solve X in primal space so X does not carry derivative information.
    to reobtain the derivatives, we evaluate michelsen's Q function instead.

    there are two parts of this function: Q1 (carries derivative information via Δ) and
    Q2 (only affects the primal value of a_assoc, not the derivatives)

    this is not necessary to do in the exact solver, as we calculate X via elementary operations that
    propagate the derivatives.
    =#
    sites = getsites(model)
    n = sites.n_sites

    Q2 = zero(first(X.v)) |> primalval
    for i ∈ @comps
        ni = n[i]
        zi = z[i]
        iszero(length(ni)) && continue
        iszero(zi) && continue
        Xᵢ = X[i]
        resᵢₐ = zero(Q2)
        for (a,nᵢₐ) ∈ pairs(ni)
            Xᵢₐ = primalval(Xᵢ[a])
            resᵢₐ += nᵢₐ * (log(Xᵢₐ) + 1 - Xᵢₐ)
        end

        Q2 += resᵢₐ*z[i]
    end
    Q1 = zero(eltype(Δ.values))
    Vinv = 1/V
    if !iszero(Vinv)
        for (idx,(i,j),(a,b)) in indices(Δ)
            Xia,nia = primalval(X[i][a]),n[i][a]
            Xjb,njb = primalval(X[j][b]),n[j][b]
            zi,zj = z[i],z[j]
            if !iszero(zi) && !iszero(zj)
                Q1 -= z[i]*z[j]*nia*njb*Xia*Xjb*(Δ.values[idx]*N_A)
            end
        end
        Q1 = Q1*Vinv
    end
    Q = Q1 + Q2
    return Q/sum(z)
end

#=
this method is used when X does propagate derivative information.

electrolyte EoS normally use this as:
_X = @f(X)
a_assoc = @f(a_assoc_impl,X)
#do something else with _X

this was the default before.
=#
function a_assoc_impl(model::EoSModel, V, T, z, X)
    _0 = zero(first(X.v))
    sites = getsites(model)
    n = sites.n_sites
    res = _0
    for i ∈ @comps
        ni = n[i]
        zi = z[i]
        iszero(zi) && continue
        iszero(length(ni)) && continue

        Xᵢ = X[i]
        resᵢₐ = _0
        for (a,nᵢₐ) ∈ pairs(ni)
            Xᵢₐ = Xᵢ[a]
            resᵢₐ +=  nᵢₐ * (log(Xᵢₐ) - Xᵢₐ*0.5 + 0.5)
        end
        res += resᵢₐ*zi
    end
    return res/sum(z)
end

#exact calculation of a_assoc when there is only one site pair
#in this case the fraction of non-bonded sites is simply xia and xjb
#so whe don't need to allocate the X vector
function a_assoc_exact_1(model::EoSModel,V,T,z,data = nothing)
    xia,xjb,i,j,a,b,n,idxs = _X_exact1(model,V,T,z,data)
    _0 = zero(xia)
    sites = getsites(model)
    nn = sites.n_sites
    res = _0
    resᵢₐ = _0
    nia = nn[i][a]
    njb = nn[j][b]
    res = z[i]*nia*(log(xia) - xia*0.5 + 0.5)
    if (i != j) | (a != b) #we check if we have 2 sites or just 1
        res += z[j]*njb*(log(xjb) - xjb*0.5 + 0.5)
    end
    return res/sum(z)
end

"""
    @assoc_loop(Xold,Xnew,expr)
Solves an association problem, given an expression for the calculation of the fraction of non-bonded sites `X`.
The macro takes care of creating the appropiate shaped vectors, and passing the appropiate iteration parameters from `AssocOptions`
Expects the following variable names in scope:
- `model` : EoS Model used
- `V`,`T`,`z` : Total volume, Temperature, mol amounts
`Xold` and `Xnew` are Vectors of Vectors, that can be indexed by component and site (`X[i][a]`).
## Example
```julia
function X(model::DAPTModel, V, T, z)
    _1 = one(V+T+first(z))
    σ = model.params.sigma.values[1][1]
    θ_c = model.params.theta_c.values[1,1][2,1]
    κ = (1 - cos(θ_c*π/180))^2/4
    ε_as = model.params.epsilon_assoc.values[1,1][2,1]
    f = exp(ε_as/(T))-1
    ρ = N_A*∑(z)/V
    Irc = @f(I)
    Xsol = @association_loop X_old X_new for i ∈ @comps, a ∈ @sites(i)
            X4 = (1-X_old[i][a])^4
            c_A = 8*π*κ*σ^3*f*(ρ*X_old[i][a]*(Irc*(1-X4) + X4/(π*ρ*σ^3)) + 2*ρ*(X_old[i][a]^2)*((1 - X_old[i][a])^3)*(Irc - 1/(π*ρ*σ^3)) )
            X_new[i][a] =1/(1+c_A)
    end
    return Xsol
end
```
"""
macro assoc_loop(Xold::Symbol,Xnew::Symbol,expr)
    return quote
        __sites = Clapeyron.getsites(model)
        idxs = __sites.n_sites.p
        X0 = fill(one(Base.promote_eltype(model,V,T,z)),length(__sites.n_sites.v))

        function x_assoc_iter!(__X_new_i,__X_old_i)
            $Xold = Clapeyron.PackedVofV(idxs,__X_old_i)
            $Xnew = Clapeyron.PackedVofV(idxs,__X_new_i)
            $expr
            return __X_new_i
        end

        options = model.assoc_options
        atol = options.atol
        rtol = options.rtol
        max_iters = options.max_iters
        α = options.dampingfactor

        Xsol = Clapeyron.Solvers.fixpoint(x_assoc_iter!,X0,Clapeyron.Solvers.SSFixPoint(α),atol=atol,rtol = rtol,max_iters = max_iters)
        Xsol
    end |> esc
end

function X_exact2!(K,X)
    k1 = K[1,2]
    k2 = K[2,1]
    #this computation is equivalent to the one done in X_exact1
    _a = k2
    _b = 1 - k2 + k1
    _c = -1
    denom = _b + sqrt(_b*_b - 4*_a*_c)
    x1 = -2*_c/denom
    x1k = k2*x1
    x2 = (1- x1k)/(1 - x1k*x1k)
    X[1] = x1
    X[2] = x2
    return X
end

function X_exact2_denseK!(K,X)
    #=
    strategy: reformulate association problem as an hyperbola
    ((x2 - t2)/b2)^2 - ((x1 - t1)/b1)^2 = 1
    x2 - t2 = t2 + b2*cosh(t)
    x1 - t1 = t1 + b1*sinh(t)

    Then we solve in exp(t) coordinates
    =#
    A1,A3,A2,A4 = K
    a = -A3/A2
    w = 1 + a + 1/(4A4) + a/(4A1)
    invert = w < 0
    
    if invert
        #if the sign of w is negative, than means that x2 has the negative sign in the hyperbola instead of x1
        A1,A4 = A4,A1
        A3,A2 = A2,A3
        a = -A3/A2
        w = 1 + a + 1/(4A4) + a/(4A1)
    end

    b1 = sqrt(-w/(a*A1))
    b2 = sqrt(w/A4)
    t1 = -0.5/A1
    t2 = -0.5/A4
    ymin = -t1/b1
    y0 = ymin
    y1 = ymin
    et = ymin + sqrt(ymin*ymin + 1)
    ε = 1e-12*one(y1)
    for i in 1:100
        x1i = t1 + b1*(et*et - 1)/(2et)
        x2i = t2 + b2*(et*et + 1)/(2et)
        y1_new = (1/(A1*x1i + A2*x2i + 1) - t1)/b1
        if y1_new < ymin
            y1_new = 0.5*ymin + 0.5*y1
        end
        y1 = y1_new
        et0 =  et
        et = 0.8*(y1 + sqrt(y1*y1 + 1)) + 0.2*et0
        abs(et - et0) < ε && break
    end
    #sinhtt = (et*et - 1)/(2et)
    #coshtt = (et*et + 1)/(2et)
    x1 = t1 + b1*(et*et - 1)/(2et)
    x2 = t2 + b2*(et*et + 1)/(2et)
    if invert
        x2,x1 = x1,x2
    end
    X[1] = x1
    X[2] = x2
    return X
end

recombine_assoc!(model) = recombine_assoc!(model,model.params.sigma)

function recombine_assoc!(model,sigma)
    _assoc_options = assoc_options(model)
    iszero(assoc_pair_length(model)) && return model
    epsilon_assoc = model.params.epsilon_assoc
    bondvol = model.params.bondvol
    bondvol,epsilon_assoc = assoc_mix(bondvol,epsilon_assoc,sigma,_assoc_options,model.sites) #combining rules for association
    copyto!(model.params.epsilon_assoc,epsilon_assoc)
    copyto!(model.params.bondvol,bondvol)
    return model
end

function assoc_matrix_solve_pure(K,idx,options)
    n = LinearAlgebra.checksquare(K) #size
    #initialization procedure:
    j = length(idx) - 1
    X = Vector{T}(undef,n)
    for i in 1:j
        v = idx[i]:idx[i+1]-1
        Xi = @view X[v]
        Ki = @view K[v,v]
        _,success = assoc_matrix_x0!(Ki,Xi)
        if !success
            res = assoc_matrix_solve(copy(Ki),options)
            Xi .= res
        end
        
    end
    return X
end
#=

=#
#=
function AX!(output,input,pack_indices,delta::Compressed4DMatrix{TT,VV} ,modelsites,ρ,z) where {TT,VV}
    _0 = zero(TT)
    p = modelsites.p::Vector{Int}
    _ii::Vector{Tuple{Int,Int}} = delta.outer_indices
    _aa::Vector{Tuple{Int,Int}} = delta.inner_indices
    _Δ::VV = delta.values
    _idx = 1:length(_ii)
    #n = modelsites
    _n::Vector{Int} = modelsites.v
    #pv.p[i]:pv.p[i+1]-1)
    @inbounds for i ∈ 1:length(z) #for i ∈ comps
        sitesᵢ = 1:(p[i+1] - p[i]) #sites are normalized, with independent indices for each component
        for a ∈ sitesᵢ #for a ∈ sites(comps(i))
            ∑X = _0
            ia = compute_index(pack_indices,i,a)
            for idx ∈ _idx #iterating for all sites
                ij = _ii[idx]
                ab = _aa[idx]
                if issite(i,a,ij,ab)
                    j = complement_index(i,ij)
                    b = complement_index(a,ab)
                    jb = compute_index(pack_indices,j,b)
                    njb = _n[jb]
                    ∑X += ρ*njb*z[j]*input[jb]*_Δ[idx]
                end
            end
            output[ia] = ∑X
        end
    end
    return output
end
=#
#res = ∑(z[i]*∑(n[i][a] * (log(X_[i][a]) - X_[i][a]/2 + 0.5) for a ∈ @sites(i)) for i ∈ @comps)/sum(z)

#=
on one site:
Xia = 1/(1+*nb*z[j]*rho*Δ*Xjb)
Xjb = 1/(1+*na*z[i]*rho*Δ*Xia)

kia = na*z[i]*rho*Δ
kjb = nb*z[j]*rho*Δ

Xia = 1/(1+kjb*Xjb)
Xjb = 1/(1+kia*Xia)

Xia = 1/(1+kjb*(1/(1+kia*Xia)))
Xia = 1/(1+kjb/(1+kia*Xia))
Xia = 1/((1+kia*Xia+kjb)/(1+kia*Xia))
Xia = (1+kia*Xia)/(1+kia*Xia+kjb)
Xia*(1+kia*Xia+kjb) = 1+kia*Xia #x = Xia
x*(1+kia*x+kjb) = 1+kia*x
x + kia*x*x + kjb*x - 1 - kia*x = 0
kia*x*x + x(kjb-kia+1) - 1 = 0
x = - (kjb-kia+1) +

x = 1/1+kiax
x(1+kx) - 1 = 0
kx2 +x - 1 = 0
end
=#

#=
function sparse_assoc_site_matrix(model,V,T,z,data=nothing)
    if data === nothing
        delta = @f(Δ)
    else
        delta = @f(Δ,data)
    end
    _sites = model.sites.n_sites
    p = _sites.p
    ρ = N_A/V
    _ii::Vector{Tuple{Int,Int}} = delta.outer_indices
    _aa::Vector{Tuple{Int,Int}} = delta.inner_indices
    _idx = 1:length(_ii)
    _Δ= delta.values
    TT = eltype(_Δ)
    count = 0
    @inbounds for i ∈ 1:length(z) #for i ∈ comps
        sitesᵢ = 1:(p[i+1] - p[i]) #sites are normalized, with independent indices for each component
        for a ∈ sitesᵢ #for a ∈ sites(comps(i))
            #ia = compute_index(pack_indices,i,a)
            for idx ∈ _idx #iterating for all sites
                ij = _ii[idx]
                ab = _aa[idx]
                issite(i,a,ij,ab) && (count += 1)
            end
        end
    end
    c1 = zeros(Int,count)
    c2 = zeros(Int,count)
    val = zeros(TT,count)
    _n = model.sites.n_sites.v
    count = 0
    @inbounds for i ∈ 1:length(z) #for i ∈ comps
        sitesᵢ = 1:(p[i+1] - p[i]) #sites are normalized, with independent indices for each component
        for a ∈ sitesᵢ #for a ∈ sites(comps(i))
            ia = compute_index(p,i,a)
            for idx ∈ _idx #iterating for all sites
                ij = _ii[idx]
                ab = _aa[idx]
                if issite(i,a,ij,ab)
                    j = complement_index(i,ij)
                    b = complement_index(a,ab)
                    jb = compute_index(p,j,b)
                    njb = _n[jb]
                    count += 1
                    c1[count] = ia
                    c2[count] = jb
                    val[count] = ρ*njb*z[j]*_Δ[idx]
                end
            end
        end
    end
    K::SparseMatrixCSC{TT,Int} = sparse(c1,c2,val)
    return K
end
#Mx = a + b(x,x)
#Axx + x - 1 = 0
#x = 1 - Axx
=#
