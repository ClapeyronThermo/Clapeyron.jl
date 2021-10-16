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
    for i ∈ 1:length(modelsites.components)
        sitesᵢ = modelsites.i_sites[i]
        iszero(length(sitesᵢ)) && continue
        for a ∈ sitesᵢ
            ∑X = _0
            @inbounds for (idx,ij,ab) in indices(delta)
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

function X(model::Union{SAFTModel,CPAModel}, V, T, z)
    _1 = one(V+T+first(z))
    ρ = N_A/V
    tol = model.absolutetolerance
    X_ = PackedVectorsOfVectors.packed_ones(typeof(_1),length(@sites(i)) for i ∈ @comps)
    idxs = indices(X_)    
    X0 = X_.v
    _Δ = @f(Δ)
    fX(out,in) = X!(out,in,idxs,_Δ,model.sites,ρ,z)
    Xsol = Solvers.fixpoint(fX,X0,Solvers.SSFixPoint(0.5),atol=tol,max_iters = 500)
    return PackedVofV(idxs,Xsol)
end

function X(model::Union{SAFTModel,CPAModel}, V, T, z,data)
    _1 = one(V+T+first(z))
    ρ = N_A/V
    tol = model.absolutetolerance
    X_ = PackedVectorsOfVectors.packed_ones(typeof(_1),length(@sites(i)) for i ∈ @comps)
    idxs = indices(X_)    
    X0 = X_.v
    _Δ = @f(Δ,data)
    fX(out,in) = X!(out,in,idxs,_Δ,model.sites,ρ,z)
    Xsol = Solvers.fixpoint(fX,X0,Solvers.SSFixPoint(0.5),atol=tol,max_iters = 500)
    return PackedVofV(idxs,Xsol)
end

function a_assoc(model::Union{SAFTModel,CPAModel}, V, T, z)
    _0 = zero(V+T+first(z))
    if iszero(length(model.sites.flattenedsites))
        return _0
    end
    X_ = @f(X)
    n = model.sites.n_sites
    res = _0
    resᵢₐ = _0
    for i ∈ @comps
        sitesᵢ = @sites(i)
        iszero(length(sitesᵢ)) && continue
        resᵢₐ = _0
        for a ∈ sitesᵢ
            resᵢₐ += n[i][a] * (log(X_[i][a]) - X_[i][a]/2 + 0.5)
        end
        res += resᵢₐ*z[i] 
    end
    return res/sum(z)
end

function a_assoc(model::Union{SAFTModel,CPAModel}, V, T, z,data)
    _0 = zero(V+T+first(z))
    if iszero(length(model.sites.flattenedsites))
        return _0
    end
    X_ = @f(X,data)
    n = model.sites.n_sites
    res = _0
    resᵢₐ = _0
    for i ∈ @comps
        sitesᵢ = @sites(i)
        iszero(length(sitesᵢ)) && continue
        resᵢₐ = _0
        for a ∈ sitesᵢ
            resᵢₐ += n[i][a] * (log(X_[i][a]) - X_[i][a]/2 + 0.5)
        end
        res += resᵢₐ*z[i] 
    end
    return res/sum(z)
end


#res =  ∑(z[i]*∑(n[i][a] * (log(X_[i][a]) - X_[i][a]/2 + 0.5) for a ∈ @sites(i)) for i ∈ @comps)/sum(z)
#=

aa = 4.76197468041569
aa = 4.76197468041569
aa = 42.788399307705504
aa = 42.788399307705504
=#