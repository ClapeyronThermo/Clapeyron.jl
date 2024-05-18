function a_assoc(model::EoSModel, V, T, z,data=nothing)
    _0 = zero(V+T+first(z))
    nn = assoc_pair_length(model)
    iszero(nn) && return _0
    isone(nn) && return a_assoc_exact_1(model,V,T,z,data)
    X,Δ = @f(X_and_Δ,data)
    return @f(a_assoc_impl,X,Δ)
end

"""
    assoc_pair_length(model::EoSModel)

Indicates the number of pair combinations between the different sites in an association model.

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
    return length(model.params.bondvol.values.values)
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

function Δ(model::EoSModel, V, T, z,data)
    Δout = assoc_similar(model,typeof(V+T+first(z)))
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

function assoc_site_matrix(model,V,T,z,ad::VV = Val{false}(),data = nothing,delta = @f(__delta_assoc,data)) where VV
    options = assoc_options(model)
    if !options.dense
        @warn "using sparse matrices for association is deprecated."
    end
    return dense_assoc_site_matrix(model,V,T,z,ad,data,delta)
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

function maybe_ad(x::X,::Val{true}) where X
    return x
end

function maybe_ad(x::X,::Val{false}) where X
    return primalval(x)
end

function dense_assoc_site_matrix(model,V,T,∂z,ad::VV = Val{false}(),data=nothing,delta = @f(__delta_assoc,data)) where VV
    sitesparam = getsites(model)
    _sites = sitesparam.n_sites
    p = _sites.p
    ρ = maybe_ad(N_A/V,ad)
    #ρ = N_A/V
    z = maybe_ad(∂z,ad)
    #z = ∂z
    _ii::Vector{Tuple{Int,Int}} = delta.outer_indices
    _aa::Vector{Tuple{Int,Int}} = delta.inner_indices
    _idx = 1:length(_ii)
    Δ = maybe_ad(delta.values,ad)
    #_Δ = delta.values
    TT = eltype(Δ)
    _n = sitesparam.n_sites.v
    nn = length(_n)
    K  = zeros(TT,nn,nn)
    count = 0
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
                    count += 1
                    K[ia,jb]  = ρ*njb*z[j]*Δ[idx]
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
    X,Δ = X_and_Δ(model,V,T,z,data,Val{true}())
    return X
end

function X_and_Δ(model::EoSModel, V, T, z,data = nothing,ad::VV = Val{false}()) where VV
    nn = assoc_pair_length(model)
    #isone(nn) && return X_exact1(model,V,T,z,data)
    options = assoc_options(model)
    _Δ = __delta_assoc(model,V,T,z,data)
    K = assoc_site_matrix(model,V,T,z,ad,data,_Δ)
    sitesparam = getsites(model)
    idxs = sitesparam.n_sites.p
    Xsol = assoc_matrix_solve(K,options)
    return PackedVofV(idxs,Xsol),_Δ
end

function assoc_matrix_solve(K,options::AssocOptions)
    atol = options.atol
    rtol = options.rtol
    max_iters = options.max_iters
    α = options.dampingfactor
    return assoc_matrix_solve(K, α, atol ,rtol, max_iters)
end

#TODO: define implicit AD here
function assoc_matrix_solve(K, α, atol ,rtol, max_iters)
    n = LinearAlgebra.checksquare(K) #size
    #initialization procedure:
    Kmin,Kmax = nonzero_extrema(K) #look for 0 < Amin < Amax
    if Kmax > 1
        f = true/Kmin
    else
        f = true-Kmin
    end
    X0 = fill(f,n) #initial point
    #=
    function to solve
    find vector x that satisfies:
    (A*x .* x) + x - 1 = 0
    Jacobian: Diagonal(Ax + 1) + Diagonal(x)*A|

    solved by reformulating in succesive substitution:
    x .= 1 ./ (1 .+ A*x)
    #we perform a "partial multiplication". that is, we use the already calculated
    #values of the next Xi to calculate the current Xi. this seems to accelerate the convergence
    #by around 50%
    =#
    function fX(out,in)
        xnorm = zero(eltype(out))
        @inbounds for i in 1:length(out)
            xnorm = min(xnorm,abs(in[i]-out[i]))
        end
        i_solved = 0
        len = length(in)
        for ii in 1:len
            Kxi = zero(eltype(out))
            @inbounds for vv in 1:i_solved
                xin = in[vv]
                xout = out[vv]
                xi = ifelse(abs(xin - xout) < xnorm, xout, xin)
                Kxi += xi*K[ii,vv]
            end
            @inbounds for vv in (i_solved+1):len
                Kxi += in[vv]*K[ii,vv]
            end
            out[ii] = 1/(1+Kxi)
            i_solved += 1
        end
        #=
        mul!(out,K,in)
        for i in 1:length(out)
            #
            Xi = in[i]
            Kxi = out[i]
            out[i] = 1/(1+Kxi)
        end 
        =#
        return out
    end

    #successive substitution until convergence
    return Solvers.fixpoint(fX,X0,Solvers.SSFixPoint(α),atol=atol,rtol = rtol,max_iters = max_iters)
end

#exact calculation of site non-bonded fraction when there is only one site

function X_exact1(model,V,T,z,data = nothing)
    xia,xjb,i,j,a,b,n,idxs = _X_exact1(model,V,T,z,data)
    pack_X_exact1(xia,xjb,i,j,a,b,n,idxs)
end

function _X_exact1(model,V,T,z,data=nothing)
    κ = model.params.bondvol.values
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
    return xia,xjb,i,j,a,b,n,idxs
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
    Q1 = zero(eltype(Δ.values))
    for (idx,(i,j),(a,b)) in indices(Δ)
        Xia,nia = primalval(X[i][a]),n[i][a]
        Xjb,njb = primalval(X[j][b]),n[j][b]
        dQ1 = z[i]*z[j]*nia*njb*Xia*Xjb*(Δ.values[idx]*N_A)
        Q1 -= dQ1
    end

    Q1 = Q1/V
    Q2 = zero(first(X.v)) |> primalval
    for i ∈ @comps
        ni = n[i]
        iszero(length(ni)) && continue
        Xᵢ = X[i]
        resᵢₐ = zero(Q2)
        for (a,nᵢₐ) ∈ pairs(ni)
            Xᵢₐ = primalval(Xᵢ[a])
            resᵢₐ += nᵢₐ * (log(Xᵢₐ) + 1 - Xᵢₐ)
        end
        Q2 += resᵢₐ*z[i]
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
        iszero(length(ni)) && continue
        Xᵢ = X[i]
        resᵢₐ = _0
        for (a,nᵢₐ) ∈ pairs(ni)
            Xᵢₐ = Xᵢ[a]
            resᵢₐ +=  nᵢₐ* (log(Xᵢₐ) - Xᵢₐ*0.5 + 0.5)
        end
        res += resᵢₐ*z[i]
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
        __sites = model.sites
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