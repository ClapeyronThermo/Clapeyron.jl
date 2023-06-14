const FIJ_TYPE = Clapeyron.PairParameter{Float64, SparseArrays.SparseMatrixCSC{Float64, Int64}}

struct GERGDepartureParam <: EoSParam
    nij::PackedSparsePairParam{Float64} #SparsePairParam #done
    tij::PackedSparsePairParam{Float64} #SparsePairParam #done
    dij::PackedSparsePairParam{Int} #SparsePairParam #done
    Fij::FIJ_TYPE #SparsePairParam #done
    beta_ij::PackedSparsePairParam{Float64} #SparsePairParam
    gamma_ij::PackedSparsePairParam{Float64} #SparsePairParam
    eta_ij::PackedSparsePairParam{Float64} #SparsePairParam
    epsilon_ij::PackedSparsePairParam{Float64} #SparsePairParam
end

@newmodelsimple GERGDeparture MultiFluidDepartureModel GERGDepartureParam

function multiparameter_a_res(model,V,T,z,departure::GERGDeparture,δ,τ,∑z = sum(z)) 
    lnδ = log(δ)
    lnτ = log(τ)
    aᵣ = multiparameter_a_res0(model,V,T,z,δ,τ,lnτ,∑z)
    _0 = zero(aᵣ)

    isone(length(z)) && return aᵣ
    F = model.pair.Fij.values
    iszero(nnz(F)) && return aᵣ
   
    Fij = nonzeros(F)
    n = model.pair.nij.values.storage
    η = model.pair.eta_ij.values.storage
    rows = rowvals(F)
    k_all = n.p
    k_exp = η.p
    nᵢⱼ = n.v
    tᵢⱼ = model.pair.tij.values.storage.v
    dᵢⱼ = model.pair.dij.values.storage.v
    ηᵢⱼ = η.v
    εᵢⱼ = model.pair.epsilon_ij.values.storage.v
    βᵢⱼ = model.pair.beta_ij.values.storage.v
    γᵢⱼ = model.pair.gamma_ij.values.storage.v
    @inbounds for j ∈ @comps
        for ii ∈ nzrange(F, j)
            i = rows[ii]
            Fᵢⱼ= Fij[ii]
            aij = _0
            #GERG2008 uses a very particular set of terms
            #instead of -η(δ-ε)^2 - β(τ-γ)^2
            #uses -η(δ-ε)^2 - β(δ-γ)
            k1,k2,kgerg = ith_index(k_all,k_exp,ii)
            
            n_pol = view(nᵢⱼ,k1)
            t_pol = view(tᵢⱼ,k1)
            d_pol = view(dᵢⱼ,k1)
            aij += term_ar_pol(δ,τ,lnδ,lnτ,_0,n_pol,t_pol,d_pol)
            
            n_gauss = view(nᵢⱼ,k2)
            t_gauss = view(tᵢⱼ,k2)
            d_gauss = view(dᵢⱼ,k2)
            η = view(ηᵢⱼ,kgerg)
            β = view(βᵢⱼ,kgerg)
            γ = view(γᵢⱼ,kgerg)
            ε = view(εᵢⱼ,kgerg)
            aij += term_ar_gerg2008(δ,τ,lnδ,lnτ,_0,n_gauss,t_gauss,d_gauss,η,β,γ,ε)
            
            aᵣ +=z[i]*z[j]*Fᵢⱼ*aij
        end
     end
    return aᵣ
end











function multiparameter_a_res(model::EmpiricMultiFluid,V,T,z,departure::GERGDeparture,δ = @f(v_scale),τ = @f(T_scale),∑z = sum(z))
    aᵣᵢ = multiparameter_a_res0(model,V,T,z,δ,τ,∑z)
    aᵣᵢⱼ = 0
    return aᵣᵢ + aᵣᵢⱼ
end