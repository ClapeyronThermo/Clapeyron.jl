function  Δ(model::Union{SAFTModel,CPAModel}, V, T, z)
    κ = model.params.bondvol.values
    Δres = zero_assoc(κ,typeof(V+T+first(z)))
    for (idx,(i,j),(a,b)) in indices(Δres)
        Δres[idx] =@f(Δ,i,j,a,b)
    end
    return Δres
end

function  Δ(model::Union{SAFTModel,CPAModel}, V, T, z,data)
    κ = model.params.bondvol.values
    Δres = zero_assoc(κ,typeof(V+T+first(z)))
    for (idx,(i,j),(a,b)) in indices(Δres)
        Δres[idx] =@f(Δ,i,j,a,b,data)
    end
    return Δres
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
    i = findfirst(>=(o-1),idxs)
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


function assoc_site_matrix(model,V,T,z,data = nothing)
    if model.assoc_options.combining == :sparse_nocombining
        K = sparse_assoc_site_matrix(model,V,T,z,data)
    else
        K = dense_assoc_site_matrix(model,V,T,z,data)
    end
    return K
end

function dense_assoc_site_matrix(model,V,T,z,data=nothing)
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
    
    _n = model.sites.n_sites.v
    
    nn = length(_n)
    K  = zeros(TT,nn,nn)
    count = 0
    combining = model.assoc_options.combining == :elliott

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
                    K[ia,jb]  = ρ*njb*z[j]*_Δ[idx]
                end
            end
        end
    end
    
    if combining
        @inbounds for ia ∈ 1:nn
            i,a = inverse_index(p,ia)
            nia = _n[ia]
            for jb ∈ 1:ia
                if iszero(K[ia,jb])               
                    j,b = inverse_index(p,jb)
                    njb = _n[jb]               
                    Δijab = sqrt(delta[i,i][a,b] * delta[j,j][a,b]) #elliott rule
                    if !iszero(Δijab)
                        K[ia,jb]  = ρ*njb*z[j]*Δijab
                        K[jb,ia]  = ρ*nia*z[i]*Δijab
                    end
                end
            end
        end
    end

    return K
end

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
    K = sparse(c1,c2,val)
    return K
end
#Mx = a + b(x,x)
#Axx + x - 1 = 0
#x = 1 - Axx

function X(model::Union{SAFTModel,CPAModel}, V, T, z,data = nothing)
    bv = model.params.bondvol.values
    nn = length(bv.values)
    isone(nn) && return X_exact1(model,V,T,z,data)
    _1 = one(V+T+first(z))
    TT = typeof(_1)    
    options = model.assoc_options
    atol = options.atol
    rtol = options.rtol
    max_iters = options.max_iters
    α = options.dampingfactor

    if model.assoc_options.combining == :sparse_nocombining
        K = sparse_assoc_site_matrix(model,V,T,z,data)
    else
        K = dense_assoc_site_matrix(model,V,T,z,data)
    end
    Kmin,Kmax = nonzero_extrema(K)
    if Kmax > 1
        f = _1/Kmin
    else
        f = _1-Kmin
    end
    
    idxs = model.sites.n_sites.p
    n = length(model.sites.n_sites.v)
    X0 = fill(f,n)

    function fX(out,in)
        mul!(out,K,in) 
        for i in 1:length(out)
            Kx = out[i]
            out[i] = _1/(_1+Kx) 
        end
        return out
    end

    Xsol::Vector{TT} = Solvers.fixpoint(fX,X0,Solvers.SSFixPoint(α),atol=atol,rtol = rtol,max_iters = max_iters)
    return PackedVofV(idxs,Xsol)
end

#exact calculation of site non-bonded fraction when there is only one site
function X_exact1(model,V,T,z,data=nothing)  
    κ = model.params.bondvol.values
    i,j = κ.outer_indices[1]
    a,b = κ.inner_indices[1]
    
    if data === nothing
        _Δ = @f(Δ,i,j,a,b)
    else
        _Δ = @f(Δ,i,j,a,b,data)
    end
    _1 = one(eltype(_Δ))
    idxs = model.sites.n_sites.p
    n = length(model.sites.n_sites.v)
    Xsol = fill(_1,n)
    _X = PackedVofV(idxs,Xsol)
    ρ = N_A/V
    zi = z[i]
    zj = z[j]
    ni = model.sites.n_sites[i]
    na = ni[a]
    nj = model.sites.n_sites[j]
    nb = nj[b]
    ρ = N_A/V
    kia = na*zi*ρ*_Δ
    kjb = nb*zj*ρ*_Δ
    #kia*x*x + x(kjb-kia+1) - 1 = 0
    _a = kia
    _b = _1 -kia + kjb
    _c = -_1
    xia = -2*_c/(_b + sqrt(_b*_b - 4*_a*_c))
    xjb = _1/(1+kia*xia)
    _X[j][b] = xjb
    _X[i][a] = xia
    return _X
end

function a_assoc(model::Union{SAFTModel,CPAModel}, V, T, z,data=nothing)
    _0 = zero(V+T+first(z))
    nn = length(model.params.bondvol.values.values)
    iszero(nn) && return _0
    X_ = @f(X,data)
    return @f(_a_assoc,X_)
end

function _a_assoc(model::Union{SAFTModel,CPAModel}, V, T, z,X_)
    _0 = zero(first(X_.v))
    n = model.sites.n_sites
    res = _0
    resᵢₐ = _0
    for i ∈ @comps
        ni = n[i]
        iszero(length(ni)) && continue        
        Xᵢ = X_[i]
        resᵢₐ = _0
        for (a,nᵢₐ) ∈ pairs(ni)
            Xᵢₐ = Xᵢ[a]
            nᵢₐ = ni[a]
            resᵢₐ +=  nᵢₐ* (log(Xᵢₐ) - Xᵢₐ/2 + 0.5)
        end
        res += resᵢₐ*z[i] 
    end
    return res/sum(z)
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
#res =  ∑(z[i]*∑(n[i][a] * (log(X_[i][a]) - X_[i][a]/2 + 0.5) for a ∈ @sites(i)) for i ∈ @comps)/sum(z)

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
