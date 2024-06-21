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

function __delta_assoc(model,V,T,z,data::M) where M
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

nonzero_extrema(K::SparseArrays.SparseMatrixCSC) = extrema(K.nzval)

function nonzero_extrema(K)
    _0 = zero(eltype(K))
    _max = _0
    _min = _0
    for k in K
        _max = max(k,_max)
        if iszero(_min)
            _min = k
        else
            if !iszero(k)
            _min = min(_min,k)
            end
        end
    end
    return _min,_max
end

function assoc_site_matrix(model,V,T,z,data = nothing,delta = @f(__delta_assoc,data))
    options = assoc_options(model)
    if !options.dense
        @warn "using sparse matrices for association is deprecated."
    end
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

function dense_assoc_site_matrix(model,V,T,z,data=nothing,delta = @f(__delta_assoc,data))
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
    _Δ = __delta_assoc(model,V,T,z,data)
    #K = assoc_site_matrix(model,primalval(V),T,primalval(z),data,primalval(_Δ))
    K = assoc_site_matrix(model,V,T,z,data,_Δ)
    sitesparam = getsites(model)
    idxs = sitesparam.n_sites.p
    Xsol = assoc_matrix_solve(K,options)
    return PackedVofV(idxs,Xsol),_Δ
end

function assoc_matrix_solve(K::AbstractMatrix{T},options::AssocOptions = AssocOptions()) where T
    atol = T(options.atol)
    rtol = T(options.rtol)
    max_iters = options.max_iters
    α = T(options.dampingfactor)
    return assoc_matrix_solve(K, α, atol ,rtol, max_iters)
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
    elseif check_antidiagonal2(K)
        X_exact2!(K,X)
        init = true
        success = true
        init = true
    elseif check_antidiagonal22(K)
    #nb-nb association with cross-association
    K11 = @view(K[1:2,1:2])
    K12 = @view(K[1:2,3:4])
    K21 = @view(K[3:4,1:2])
    K22 = @view(K[3:4,3:4])
    if (iszero(K12) & iszero(K21)) | iszero(K11) | iszero(K22)
        #solve each association separately, if one of the diagonal association
        #submatrices is zero, then cross-association does not have any sense.
        X_exact2!(K11,@view(X[1:2]))
        X_exact2!(K22,@view(X[3:4]))
    else
        #general solution.
        X_exact4!(K,X)
    end
    success = true
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

function assoc_matrix_solve(K::AbstractMatrix{T}, α::T, atol ,rtol, max_iters) where T
    n = LinearAlgebra.checksquare(K) #size
    #initialization procedure:
    X0 = Vector{T}(undef,n)
    X0,success = assoc_matrix_x0!(K,X0)
    success && return X0
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
            if finite
                return Xsol
            else
                Xsol .= NaN
                return Xsol
            end
        end
        X0 .= Xsol
       # @show Xsol
    end
    H = Matrix{T}(undef,n,n)
    H .= 0
    piv = zeros(Int,n)
    F = Solvers.unsafe_LU!(H,piv)
    if !converged #proceed to newton minimization
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
            F = Solvers.unsafe_LU!(F)
            dX .= 1 ./ Xsol .- 1 .- KX #gradient
            ldiv!(F,dX) #we solve H/g, overwriting g
            for k in 1:length(dX)
                Xk = Xsol[k]
                dXk = dX[k]
                X_newton = Xk - dXk
                if !(0 <= X_newton <= 1)
                    Xsol[k] = 1/(1 + KX[k]) #successive substitution step
                else
                    Xsol[k] = X_newton #newton step
                end
            end
           # Xsol .-= dX

            converged,finite = Solvers.convergence(Xsol,X0,atol,rtol,false,Inf)
            #@show converged,finite
            if converged
                if !finite
                    fill!(Xsol,NaN)
                end
                return Xsol
            end
            X0 .= Xsol
        end
    end
    if !converged
        Xsol .= NaN
    end
    return Xsol
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
macro assoc_loop(Xold,Xnew,expr)
    return quote
        __sites = getsites(model)
        idxs = __sites.n_sites.p
        X0 = fill(one(V+T+first(z)),length(__sites.n_sites.v))

        function x_assoc_iter!(__X_new_i,__X_old_i)
            $Xold = PackedVofV(idxs,__X_old_i)
            $Xnew = PackedVofV(idxs,__X_old_i)
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

function X_exact4!(K,X)
    #=
    strategy is the following:
    given K (4x4 assoc matrix) and X

    1. we solve for x1 (7th order polynomial, bracketed newton)
    2. we solve for x3 (2nd order polynomial, analytic)
    3. x2 and x4 only depend on x1 and x3.

    this supposes non-zero values of the diagonal submatrices.
    we catch the zero case earlier.
    =#
    _0 = zero(eltype(K))
    k3 = K[2]
    k7 = K[4]
    k1 = K[5]
    k5 = K[7]
    k4 = K[10]
    k8 = K[12]
    k2 = K[13]
    k6 = K[15]

    X_exact2!(@view(K[1:2,1:2]),@view(X[1:2]))
    x10 = X[1]
    pol_x1 = __assoc_x1_poly(K)
    dpol_x1 = Solvers.polyder(pol_x1)
    function f0(x)
        fx = evalpoly(x,pol_x1)
        dfx = evalpoly(x,dpol_x1)
        return fx,fx/dfx
    end
    prob_x1 = Roots.ZeroProblem(f0,(zero(x0),one(x0),x10))
    x1 = Roots.solve(prob_x1,Roots.LithBoonkkampIJzermanBracket()) #bracketed newton
    x3 = __assoc_x3(K,x1)
    x2 = 1 / (1 + k3*x1 + k4*x3)
    x4 =  1 / (1 + k7*x1 + k8*x3)
    X[1] = x1
    X[2] = x2
    X[3] = x3
    X[4] = x4
    return X
end

function __assoc_x1_poly(K)
    k3 = K[2]
    k7 = K[4]
    k1 = K[5]
    k5 = K[7]
    k4 = K[10]
    k8 = K[12]
    k2 = K[13]
    k6 = K[15]
    var1 = k8*k8
    var2 = k1*k2*k4*k5*k8
    var3 = k1*k1
    var4 = var3*k4*k6*k8
    var5 = k4*k4
    var6 = var3*var5*k6*k8
    var7 = -(var3*k4*k5*k6*k8)
    var8 = -(2*k1*k2*k4*k5*k6*k8)
    var9 = k6*k6
    var10 = var3*k4*var9*k8
    var11 = -(k1*k2*k5*var1)
    var12 = k1*k2*k4*k5*var1
    var13 = k5*k5
    var14 = -(k1*k2*var13*var1)
    var15 = -(var3*k6*var1)
    var16 = -(var3*k4*k6*var1)
    var17 = 2*var3*k5*k6*var1
    var18 = k1*k2*k5*k6*var1
    var19 = var1*k8
    var20 = -(k1*k2*k5*var19)
    var21 = k2*k2
    var22 = k1*k2*k3*k4*k5*k8
    var23 = var3*k1
    var24 = var3*k3*k4*k6*k8
    var25 = -(2*k1*k2*k3*k4*k5*k6*k8)
    var26 = var3*k3*k4*var9*k8
    var27 = 2*k1*k2*k4*k5*k7*k8
    var28 = 2*var3*k4*k6*k7*k8
    var29 = 2*var3*var5*k6*k7*k8
    var30 = -(2*var3*k4*k5*k6*k7*k8)
    var31 = -(2*k1*k2*k4*k5*k6*k7*k8)
    var32 = var3*k4*var9*k7*k8
    var33 = k7*k7
    var34 = -(2*k1*k2*k3*k5*var1)
    var35 = k1*k2*k3*k4*k5*var1
    var36 = -(k1*k2*k3*var13*var1)
    var37 = -(2*var3*k3*k6*var1)
    var38 = -(var3*k3*k4*k6*var1)
    var39 = 2*var3*k3*k5*k6*var1
    var40 = 2*k1*k2*k3*k5*k6*var1
    var41 = k3*k3
    var42 = -(k1*k2*k5*k7*var1)
    var43 = k1*k2*k4*k5*k7*var1
    var44 = -(k1*k2*var13*k7*var1)
    var45 = -(var3*k6*k7*var1)
    var46 = -(var3*k4*k6*k7*var1)
    var47 = 2*var3*k5*k6*k7*var1
    var48 = -(2*k1*k2*k3*k5*var19)
    var49 = k2*var21
    var50 = 2*k1*k2*k3*k4*k5*k7*k8
    var51 = 2*var3*k3*k4*k6*k7*k8
    var52 = -(2*k1*k2*k3*k4*k5*k6*k7*k8)
    var53 = var3*k3*k4*var9*k7*k8
    var54 = k1*k2*k4*k5*var33*k8
    var55 = var3*k4*k6*var33*k8
    var56 = var3*var5*k6*var33*k8
    var57 = -(var3*k4*k5*k6*var33*k8)
    var58 = -(k1*k2*var41*k5*var1)
    var59 = var3*var3
    var60 = -(var3*var41*k6*var1)
    var61 = k1*k2*var41*k5*k6*var1
    var62 = -(2*k1*k2*k3*k5*k7*var1)
    var63 = k1*k2*k3*k4*k5*k7*var1
    var64 = -(k1*k2*k3*var13*k7*var1)
    var65 = -(2*var3*k3*k6*k7*var1)
    var66 = -(var3*k3*k4*k6*k7*var1)
    var67 = 2*var3*k3*k5*k6*k7*var1
    var68 = -(k1*k2*var41*k5*var19)
    var69 = k1*k2*k3*k4*k5*var33*k8
    var70 = var3*k3*k4*k6*var33*k8
    var71 = -(k1*k2*var41*k5*k7*var1)
    var72 = -(var3*var41*k6*k7*var1)
    p0 = zero(eltype(K))
    p1 = k1*k4*k5*k6*k8-k1*k5*k6*var1
    p2 =var20-k1*k5*k6*k7*var1-2*k1*k3*k5*k6*var1+var18+var17+3*k1*k5*k6*var1+var16+var15+var14+var12+var11+2*k1*k4*k5*k6*k7*k8+var10+k1*k3*k4*k5*k6*k8+var8+var7-3*k1*k4*k5*k6*k8+var6+var4+var2
    p3 = var48+var3*k2*k5*var19+2*k1*k2*k5*var19-var3*k2*k4*var19-var3*k2*var19-2*k1*k3*k5*k6*k7*var1+var47+3*k1*k5*k6*k7*var1+var46+var45+var44+var43+var42-k1*var41*k5*k6*var1+var40+var39+6*k1*k3*k5*k6*var1-var3*k2*k5*k6*var1-2*k1*k2*k5*k6*var1-var23*k5*k6*var1-4*var3*k5*k6*var1-3*k1*k5*k6*var1+var38+2*var3*k2*k4*k6*var1+var23*k4*k6*var1+2*var3*k4*k6*var1+var37+var3*k2*k6*var1+2*var23*k6*var1+2*var3*k6*var1+var36+k1*var21*var13*var1+var3*k2*var13*var1+2*k1*k2*var13*var1+var35-k1*var21*k4*k5*var1-2*var3*k2*k4*k5*var1-2*k1*k2*k4*k5*var1+var34+k1*var21*k5*var1+2*k1*k2*k5*var1+var3*k2*var5*var1-var3*k2*var1+k1*k4*k5*k6*var33*k8+var32+2*k1*k3*k4*k5*k6*k7*k8+var31+var30-6*k1*k4*k5*k6*k7*k8+var29+var28+var27+var26-var3*k2*k4*var9*k8-var23*k4*var9*k8-2*var3*k4*var9*k8+var25-3*k1*k3*k4*k5*k6*k8+k1*var21*k4*k5*k6*k8+var3*k2*k4*k5*k6*k8+4*k1*k2*k4*k5*k6*k8+2*var3*k4*k5*k6*k8+3*k1*k4*k5*k6*k8-var3*k2*var5*k6*k8-2*var3*var5*k6*k8+var24-var23*k4*k6*k8-2*var3*k4*k6*k8+var22-2*k1*var21*k4*k5*k8-var3*k2*k4*k5*k8-2*k1*k2*k4*k5*k8+var3*k2*var5*k8+var3*k2*k4*k8
    p4 = var68+var3*k2*k3*k5*var19+4*k1*k2*k3*k5*var19-var3*k2*k5*var19+var20-var3*k2*k3*k4*var19+var3*k2*k4*var19-2*var3*k2*k3*var19+var23*k2*var19+var3*k2*var19-k1*var41*k5*k6*k7*var1+var67+6*k1*k3*k5*k6*k7*var1-var23*k5*k6*k7*var1-4*var3*k5*k6*k7*var1-3*k1*k5*k6*k7*var1+var66+var23*k4*k6*k7*var1+2*var3*k4*k6*k7*var1+var65+2*var23*k6*k7*var1+2*var3*k6*k7*var1+var64+var3*k2*var13*k7*var1+2*k1*k2*var13*k7*var1+var63-2*var3*k2*k4*k5*k7*var1-2*k1*k2*k4*k5*k7*var1+var62+2*k1*k2*k5*k7*var1+var3*k2*var5*k7*var1-var3*k2*k7*var1+var61+3*k1*var41*k5*k6*var1-var3*k2*k3*k5*k6*var1-4*k1*k2*k3*k5*k6*var1-4*var3*k3*k5*k6*var1-6*k1*k3*k5*k6*var1+var3*k2*k5*k6*var1+var18+var23*k5*k6*var1+var17+k1*k5*k6*var1+2*var3*k2*k3*k4*k6*var1+2*var3*k3*k4*k6*var1-2*var3*k2*k4*k6*var1-var23*k4*k6*var1+var16+var60+2*var3*k2*k3*k6*var1+2*var23*k3*k6*var1+4*var3*k3*k6*var1-var23*k2*k6*var1-var3*k2*k6*var1-var59*k6*var1-2*var23*k6*var1+var15+k1*var21*k3*var13*var1+2*k1*k2*k3*var13*var1-k1*var21*var13*var1-var3*k2*var13*var1+var14-k1*var21*k3*k4*k5*var1-2*k1*k2*k3*k4*k5*var1+k1*var21*k4*k5*var1+2*var3*k2*k4*k5*var1+var12+var58+2*k1*var21*k3*k5*var1+4*k1*k2*k3*k5*var1+var3*var21*k5*var1-k1*var21*k5*var1+var23*k2*k5*var1+var11-var3*k2*var5*var1+var3*var21*k4*var1-var23*k2*k4*var1-2*var3*k2*k3*var1+var3*var21*var1+var23*k2*var1+var3*k2*var1+k1*k3*k4*k5*k6*var33*k8+var57-3*k1*k4*k5*k6*var33*k8+var56+var55+var54+var53-var23*k4*var9*k7*k8-2*var3*k4*var9*k7*k8+var52-6*k1*k3*k4*k5*k6*k7*k8+var3*k2*k4*k5*k6*k7*k8+4*k1*k2*k4*k5*k6*k7*k8+4*var3*k4*k5*k6*k7*k8+6*k1*k4*k5*k6*k7*k8-var3*k2*var5*k6*k7*k8-4*var3*var5*k6*k7*k8+var51-2*var23*k4*k6*k7*k8-4*var3*k4*k6*k7*k8+var50-2*k1*var21*k4*k5*k7*k8-2*var3*k2*k4*k5*k7*k8-4*k1*k2*k4*k5*k7*k8+2*var3*k2*var5*k7*k8+2*var3*k2*k4*k7*k8-var3*k2*k3*k4*var9*k8-2*var3*k3*k4*var9*k8+var3*k2*k4*var9*k8+var23*k4*var9*k8+var10+k1*var21*k3*k4*k5*k6*k8+4*k1*k2*k3*k4*k5*k6*k8+3*k1*k3*k4*k5*k6*k8-k1*var21*k4*k5*k6*k8-var3*k2*k4*k5*k6*k8+var8+var7-k1*k4*k5*k6*k8+var3*k2*var5*k6*k8+var6-2*var3*k3*k4*k6*k8-var3*var21*k4*k6*k8-var23*k2*k4*k6*k8+var23*k4*k6*k8+var4-2*k1*var21*k3*k4*k5*k8-2*k1*k2*k3*k4*k5*k8+k1*var49*k4*k5*k8+var3*var21*k4*k5*k8+2*k1*var21*k4*k5*k8+var3*k2*k4*k5*k8+var2-var3*var21*var5*k8-var3*k2*var5*k8+var3*k2*k3*k4*k8-var3*var21*k4*k8-var23*k2*k4*k8-var3*k2*k4*k8
    p5 = 2*k1*k2*var41*k5*var19-var3*k2*k3*k5*var19+var48+var3*k2*k3*k4*var19-var3*k2*var41*var19+var23*k2*k3*var19+2*var3*k2*k3*var19+3*k1*var41*k5*k6*k7*var1-4*var3*k3*k5*k6*k7*var1-6*k1*k3*k5*k6*k7*var1+var23*k5*k6*k7*var1+var47+k1*k5*k6*k7*var1+2*var3*k3*k4*k6*k7*var1-var23*k4*k6*k7*var1+var46+var72+2*var23*k3*k6*k7*var1+4*var3*k3*k6*k7*var1-var59*k6*k7*var1-2*var23*k6*k7*var1+var45+2*k1*k2*k3*var13*k7*var1-var3*k2*var13*k7*var1+var44-2*k1*k2*k3*k4*k5*k7*var1+2*var3*k2*k4*k5*k7*var1+var43+var71+4*k1*k2*k3*k5*k7*var1+var23*k2*k5*k7*var1+var42-var3*k2*var5*k7*var1-var23*k2*k4*k7*var1-2*var3*k2*k3*k7*var1+var23*k2*k7*var1+var3*k2*k7*var1-2*k1*k2*var41*k5*k6*var1-3*k1*var41*k5*k6*var1+var3*k2*k3*k5*k6*var1+var40+var39+2*k1*k3*k5*k6*var1-2*var3*k2*k3*k4*k6*var1+var38+var3*k2*var41*k6*var1+2*var3*var41*k6*var1-var23*k2*k3*k6*var1-2*var3*k2*k3*k6*var1-2*var23*k3*k6*var1+var37-k1*var21*k3*var13*var1+var36+k1*var21*k3*k4*k5*var1+var35+k1*var21*var41*k5*var1+2*k1*k2*var41*k5*var1+var3*var21*k3*k5*var1-2*k1*var21*k3*k5*var1+var34+var3*var21*k3*k4*var1-var3*k2*var41*var1+2*var3*var21*k3*var1+var23*k2*k3*var1+2*var3*k2*k3*var1-3*k1*k3*k4*k5*k6*var33*k8+2*var3*k4*k5*k6*var33*k8+3*k1*k4*k5*k6*var33*k8-2*var3*var5*k6*var33*k8+var70-var23*k4*k6*var33*k8-2*var3*k4*k6*var33*k8+var69-var3*k2*k4*k5*var33*k8-2*k1*k2*k4*k5*var33*k8+var3*k2*var5*var33*k8+var3*k2*k4*var33*k8-2*var3*k3*k4*var9*k7*k8+var23*k4*var9*k7*k8+var32+4*k1*k2*k3*k4*k5*k6*k7*k8+6*k1*k3*k4*k5*k6*k7*k8-var3*k2*k4*k5*k6*k7*k8+var31+var30-2*k1*k4*k5*k6*k7*k8+var3*k2*var5*k6*k7*k8+var29-4*var3*k3*k4*k6*k7*k8-var23*k2*k4*k6*k7*k8+2*var23*k4*k6*k7*k8+var28-2*k1*var21*k3*k4*k5*k7*k8-4*k1*k2*k3*k4*k5*k7*k8+var3*var21*k4*k5*k7*k8+2*k1*var21*k4*k5*k7*k8+2*var3*k2*k4*k5*k7*k8+var27-var3*var21*var5*k7*k8-2*var3*k2*var5*k7*k8+2*var3*k2*k3*k4*k7*k8-var3*var21*k4*k7*k8-2*var23*k2*k4*k7*k8-2*var3*k2*k4*k7*k8+var3*k2*k3*k4*var9*k8+var26-k1*var21*k3*k4*k5*k6*k8+var25-k1*k3*k4*k5*k6*k8-var3*var21*k3*k4*k6*k8+var24+k1*var49*k3*k4*k5*k8+2*k1*var21*k3*k4*k5*k8+var22-var3*var21*k3*k4*k8-var3*k2*k3*k4*k8
    p6 = var68+var3*k2*var41*var19-3*k1*var41*k5*k6*k7*var1+var67+2*k1*k3*k5*k6*k7*var1+var66+2*var3*var41*k6*k7*var1-2*var23*k3*k6*k7*var1+var65+var64+var63+2*k1*k2*var41*k5*k7*var1+var62-var3*k2*var41*k7*var1+var23*k2*k3*k7*var1+2*var3*k2*k3*k7*var1+var61+k1*var41*k5*k6*var1-var3*k2*var41*k6*var1+var60-k1*var21*var41*k5*var1+var58+var3*var21*var41*var1+var3*k2*var41*var1+3*k1*k3*k4*k5*k6*var33*k8+var57-k1*k4*k5*k6*var33*k8+var56-2*var3*k3*k4*k6*var33*k8+var23*k4*k6*var33*k8+var55-2*k1*k2*k3*k4*k5*var33*k8+var3*k2*k4*k5*var33*k8+var54-var3*k2*var5*var33*k8+var3*k2*k3*k4*var33*k8-var23*k2*k4*var33*k8-var3*k2*k4*var33*k8+var53+var52-2*k1*k3*k4*k5*k6*k7*k8+var51+2*k1*var21*k3*k4*k5*k7*k8+var50-var3*var21*k3*k4*k7*k8-2*var3*k2*k3*k4*k7*k8
    p7 = k1*var41*k5*k6*k7*var1+var72+var71+var3*k2*var41*k7*var1-k1*k3*k4*k5*k6*var33*k8+var70+var69-var3*k2*k3*k4*var33*k8
    if iszero(p7)
        return (p1,p2,p3,p4,p5,p6,p7)
    else
        return (p1/p7,p2/p7,p3/p7,p4/p7,p5/p7,p6/p7,one(eltype(K)))
    end
end

function __assoc_x3(K,x1)
    _0 = zero(eltype(K))
    k3 = K[2]
    k7 = K[4]
    k1 = K[5]
    k5 = K[7]
    k4 = K[10]
    k8 = K[12]
    k2 = K[13]
    k6 = K[15]
    c3 = -k3*k7
    c2 = k7*(k3 - k1 - 1) - k3*(k2 + 1)
    c1 = k7 + k3 - k2 - k1 - 1
    c = evalpoly(x1,(one(c1),c1,c2,c3))
    b2 = -k3*k8-k4*k7
    b1 = k8*(k3 - k1 - 1) + k4*(k7 -k2 - 1)
    b0 = k8+k4
    b = evalpoly(x1,(b0,b1,b2))
    a = k4*k8*(1 - x1)
    Δ = b*b - 4*a*c
    return (-b + sqrt(Δ))/(2*a)
end

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