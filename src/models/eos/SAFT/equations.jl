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
function issite(i,a,ij,ab)
    ia = (i,a)
    i1,i2 = ij
    a1,a2 = ab
    ia1 = (i1,a1)
    ia2 = (i2,a2)
    return (ia == ia1) | (ia == ia2)
end

function X!(output,input,pack_indices,delta,modelsites,ρ,z)
    
    _0 = zero(eltype(output))
    _1 = one(_0)
    Xinput = PackedVofV(pack_indices,input)
    Xoutput = PackedVofV(pack_indices,output)
    n = modelsites.n_sites
    for i ∈ 1:length(modelsites.components) #for i ∈ comps
        sitesᵢ = 1:length(n[i]) #sites are normalized, with independent indices for each component
        iszero(length(sitesᵢ)) && continue
        for a ∈ sitesᵢ #for a ∈ sites(comps(i))
            ∑X = _0
            for (idx,ij,ab) in indices(delta) #iterating for all sites
                if issite(i,a,ij,ab)
                    j = complement_index(i,ij)
                    b = complement_index(a,ab)
                    nj = n[j]
                    njb = nj[b]
                    ∑X += ρ*njb*z[j]*Xinput[j][b]*delta.values[idx]
                end
            end
            rhs = _1/(_1 +∑X)
            Xoutput[i][a] = rhs
        end
    end
    return output
end

function X(model::Union{SAFTModel,CPAModel}, V, T, z,data = nothing)
    _1 = one(V+T+first(z))
    X_ = PackedVectorsOfVectors.packed_ones(typeof(_1),length(@sites(i)) for i ∈ @comps)
    idxs = indices(X_)    
    X0 = X_.v
    bv = model.params.bondvol.values
    nn = length(bv.values)
    if isone(nn)
        xia,xjb = X_exact1(model,V,T,z,data)
        i,j = only(bv.outer_indices) 
        a,b = only(bv.inner_indices) 
        #in the case that i = j, a = b, this does assignment twice
        #we do xia last because xjb is calculated from xia
        X_[j][b] = xjb
        X_[i][a] = xia
        return X_
    end
    ρ = N_A/V
    if data === nothing
        _Δ = @f(Δ)
    else
        _Δ = @f(Δ,data)
    end  
    fX(out,in) = X!(out,in,idxs,_Δ,model.sites,ρ,z)

    options = model.assoc_options
    atol = options.atol
    rtol = options.rtol
    max_iters = options.max_iters
    α = options.dampingfactor
    Xsol = Solvers.fixpoint(fX,X0,Solvers.SSFixPoint(α),atol=atol,rtol = rtol,max_iters = max_iters)
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
    _b = 1 -kia + kjb
    _c = -1
    xia = -2*_c/(_b + sqrt(_b*_b - 4*_a*_c))
    xjb = 1/(1+kia*xia)
    return (xia,xjb)
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

function show_param(param)
    for i in fieldnames(typeof(param))
        display(getfield(param,i))
    end
end

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