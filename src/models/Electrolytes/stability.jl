function tpd_plan(model::ISElectrolyteModel,z,is_liquidz,lle,id_test,K_test,pure_test)
    plan = Tuple{Symbol,Symbol,NTuple{3,Int}}[]
    nc = length(model)
    neutral = ones(Bool,length(model))
    isalts = model.salt.isalts
    neutral[isalts] .= false

    if is_liquidz && id_test && !lle && iszero(length(isalts))
        push!(plan,(:ideal_gas,:vapour,0))
    end

    ids = sortperm(z)
    if is_liquidz
        for i in 1:nc
            #idx_solvent = ids[i]
            #if neutral[idx_solvent]
            #    push!(plan,(:pure,:liquid,(idx_solvent,0,0)))
            #end
            if i != nc
                push!(plan,(:pereira,:liquid,(i,-1,0)))
                push!(plan,(:pereira,:liquid,(i,1,0)))
            end
        end
        if !lle
            for i in 1:nc
                !iszero(ZZ[i]) && push!(plan,(:pure,:vapour,(ids[i],0,0)))
            end
        end
    else
        for i in 1:nc
            #idx_solvent = ids[i]
            #if neutral[idx_solvent]          
            #    push!(plan,(:pure,:liquid,(idx_solvent,0,0)))
            #end

            if i != nc
                push!(plan,(:pereira,:liquid,(i,-1,0)))
                push!(plan,(:pereira,:liquid,(i,1,0)))
            end
        end
    end
    return plan
end

function tpd_ss_ψ(ψ,K,Z)
    res = zero(Base.promote_eltype(ψ,K,Z))
    for i in 1:length(K)
      wi = K[i]*exp(Z[i]*ψ)
      res += Z[i]*wi 
    end
    
    #=
    # log(w) - log(z) = lnphi(z)/lnphi(w) + Z*ψ
    # log(w) + logphi(w) - logphi(z) - log(z) - Z*ψ = 0
    # 
    # =#
    return res
end

function __tpd_ss_update!(w,model::ESElectrolyteModel,d,z,lnϕw)
    #= dz = logphiz + log(z)
    w = exp(logphiz - logphiw + log(z))
    w = z*log(K)
    =#
    w .= exp.(d .- lnϕw) #K is stored in w
    Z = model.charge
    f(x) = tpd_ss_ψ(x,w,Z)
    prob = Roots.ZeroProblem(f,zero(Base.promote_eltype(w,d)))
    ψ = Roots.solve(prob)
    S = zero(eltype(w))
    for i in eachindex(w)
        w[i] *= exp(ψ * Z[i])
        S += w[i]
    end
    
    
    zz = sum(z)
    tpd = oneunit(S)
    K_norm = zero(S)
    for i in eachindex(w)
        wi = w[i]
        K_norm += log(zz*wi/z[i])^2
        tpd += wi*(log(wi) + lnϕw[i] - d[i] - 1)
    end
    w ./= S
    return S,tpd,K_norm
end

function tpd_plan(model::ESElectrolyteModel,z,is_liquidz,lle,id_test,K_test,pure_test)
    plan = Tuple{Symbol,Symbol,NTuple{3,Int}}[]
    nc = length(model)
    neutral = ones(Bool,length(model))
    Z = model.charge
    neutral = iszero.(Z)

    if is_liquidz && id_test && !lle && iszero(length(Z))
        push!(plan,(:ideal_gas,:vapour,(0,0,0)))
    end

    ids = sortperm(z)
    if is_liquidz
        for i in 1:nc
            #idx_solvent = ids[i]
            #if neutral[idx_solvent]
            #    push!(plan,(:pure,:liquid,(idx_solvent,0,0)))
            #end
            if i != nc
                push!(plan,(:pereira,:liquid,(i,-1,0)))
                push!(plan,(:pereira,:liquid,(i,1,0)))
            end
        end
        if !lle
            for i in 1:nc
                !iszero(ZZ[i]) && push!(plan,(:pure,:vapour,(ids[i],0,0)))
            end
        end
    else
        for i in 1:nc
            #idx_solvent = ids[i]
            #if neutral[idx_solvent]          
            #    push!(plan,(:pure,:liquid,(idx_solvent,0,0)))
            #end

            if i != nc
                push!(plan,(:pereira,:liquid,(i,-1,0)))
                push!(plan,(:pereira,:liquid,(i,1,0)))
            end
        end
    end
    return plan
end
