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
    #=
    n = modelsites.n_flattenedsites
    sitesⱼ = modelsites.i_flattenedsites
    for i ∈ 1:length(modelsites.components)
        sitesᵢ = modelsites.i_sites[i]
        iszero(length(sitesᵢ)) && continue
        idxi = sitesⱼ[i]
        for a ∈ sitesᵢ
            ∑X = _0
            for (idx,ij,ab) in indices(delta)
                if issite(i,a,ij,ab)
                    j = complement_index(i,ij)
                    b = complement_index(a,ab)
                    idxj = sitesⱼ[j]
                    _b =  idxj[b]
                    nj = n[j]
                    njb = nj[b]
                    ∑X += ρ*njb*z[j]*Xoutput[j][_b]*delta.values[idx]   
                end
            end
            _a = idxi[a]
            rhs = _1/(_1 +∑X)
            Xoutput[i][_a] = rhs
        end
    end=#
        n = modelsites.n_sites
    for i ∈ 1:length(modelsites.components)
        sitesᵢ = modelsites.i_sites[i]
        iszero(length(sitesᵢ)) && continue
        for a ∈ sitesᵢ
            ∑X = _0
            for (idx,ij,ab) in indices(delta)
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
    tol = model.absolutetolerance
    X_ = PackedVectorsOfVectors.packed_ones(typeof(_1),length(@sites(i)) for i ∈ @comps)
    idxs = indices(X_)    
    X0 = X_.v
    nn = length(model.params.bondvol.values.values)
    if isone(nn)
        x1,x2 = X_exact1(model,V,T,z,data)
        X0[1] = x1
        X0[2] = x2
        return PackedVofV(idxs,X0)
    end
    ρ = N_A/V
    if data === nothing
        _Δ = @f(Δ)
    else
        _Δ = @f(Δ,data)
    end  
    fX(out,in) = X!(out,in,idxs,_Δ,model.sites,ρ,z)
    Xsol = Solvers.fixpoint(fX,X0,Solvers.SSFixPoint(0.5),atol=tol,max_iters = 500)
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
    ni = model.sites.n_flattenedsites[i]
    na = ni[a]
    nj = model.sites.n_flattenedsites[j]
    nb = nj[b]
    ρ = N_A/V
    kia = na*zi*ρ*_Δ
    kjb = nb*zj*ρ*_Δ
    #kia*x*x + x(kjb-kia+1) - 1 = 0
    a = kia
    b = 1 -kia + kjb
    c = -1
    xia = -2*c/(b + sqrt(b*b - 4*a*c))
    xjb = 1/(1+kia*xia)
    if (i,a) < (j,b)
        return (xia,xjb)
    else
        return (xjb,xia)
    end 
end

function a_assoc(model::Union{SAFTModel,CPAModel}, V, T, z)
    _0 = zero(V+T+first(z))
    nn = length(model.params.bondvol.values.values)
    iszero(nn) && return _0
    X_ = @f(X)
    return @f(_a_assoc,X_)
end

function a_assoc(model::Union{SAFTModel,CPAModel}, V, T, z,data)
    _0 = zero(V+T+first(z))
    nn = length(model.params.bondvol.values.values)
    iszero(nn) && return _0
    X_ = @f(X,data)
    return @f(_a_assoc,X_)
end

function _a_assoc(model::Union{SAFTModel,CPAModel}, V, T, z,X_)
    _0 = zero(first(X_.v))
    n = model.sites.n_flattenedsites
    _idx = model.sites.i_flattenedsites
    res = _0
    resᵢₐ = _0
    for i ∈ @comps
        sitesᵢ = @sites(i)
        _i = _idx[i] 
        iszero(length(sitesᵢ)) && continue
        ni = n[i]
        Xᵢ = X_[i]
        resᵢₐ = _0
        for a ∈ sitesᵢ
            _a = _i[a]
            Xᵢₐ = Xᵢ[_a]
            nᵢₐ = ni[a]
            resᵢₐ +=  nᵢₐ* (log(Xᵢₐ) - Xᵢₐ/2 + 0.5)
        end
        res += resᵢₐ*z[i] 
    end
    return res/sum(z)

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

=#