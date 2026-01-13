function tpd_plan(model::ElectrolyteModel,z,is_liquidz,lle,id_test,K_test,pure_test)
    plan = Tuple{Symbol,Symbol,NTuple{3,Int}}[]

    if is_liquidz && id_test && !lle
        push!(plan,(:ideal_gas,:vapour,0))
    end

    #=
    if K_test
        if is_liquidz
            lle || push!(plan,(:K,:vapour,1))
            push!(plan,(:K,:liquid,1))
            push!(plan,(:K,:liquid,-1))
            lle || push!(plan,(:K,:vapour,-1))
        else
            push!(plan,(:K,:liquid,1))
            push!(plan,(:K,:liquid,-1))
        end
    end =#
    Z = model.charge
    _,i_elec_pivot = findmax(Z)
    if pure_test
        ids = sortperm(z)
        ZZ = @view Z[ids]
        if is_liquidz
            for i in 1:length(z)
                if iszero(ZZ[i])
                    idx_solvent = ids[i]
                    push!(plan,(:pure,:liquid,(idx_solvent,0,0)))
                    for j in 1:length(z)
                        if !iszero(ZZ[j])
                            idx_elec = ids[j]
                            xx = (idx_solvent,idx_elec,i_elec_pivot)
                            i_elec_pivot != idx_elec && push!(plan,(:electrolyte_balanced,:liquid,xx))

                        end
                    end
                end
            end
            if !lle
                for i in 1:length(z)
                    !iszero(ZZ[i]) && push!(plan,(:pure,:vapour,ids[i]))
                end
            end
        else
            for i in 1:length(z)
                if iszero(ZZ[i])
                    idx_solvent = ids[i]
                    push!(plan,(:pure,:liquid,(idx_solvent,0,0)))
                    for j in 1:length(z)
                        if !iszero(ZZ[j])
                            idx_elec = ids[j]
                            xx = (idx_solvent,idx_elec,i_elec_pivot)
                            i_elec_pivot != idx_elec && push!(plan,(:electrolyte_balanced,:liquid,xx))
                        end
                    end
                end
            end
        end
    end
    return plan
end