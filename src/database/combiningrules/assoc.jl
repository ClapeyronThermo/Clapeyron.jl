function assoc_extend(param::AssocParam{T}) where T
    vals = param.values
    length(vals.values) == 0 && return param
    extended_vals = Compressed4DMatrices.c4d_from_site_offsets(T,vals.site_offsets)
    for (idx,(i,j),(a,b)) in indices(vals)
        extended_vals[(i,j,a,b)] = vals[(i,j,a,b)]
    end
    return AssocParam(param.name,param.components,extended_vals,param.sites,param.sourcecsvs,param.sources)
end

bondvol_mix(bondvol::AssocParam) = bondvol_mix(bondvol,nothing)

function bondvol_mix(bondvol::AssocParam,::Nothing)
    length(bondvol.values.values) == 0 && return deepcopy(bondvol)

    param = assoc_extend(bondvol)
    mat = param.values
    for (idx,(i,j),(a,b)) in indices(mat)
        dij = sqrt(mat[i,i][a,b]*mat[j,j][a,b])
        if !iszero(dij) && iszero(mat[idx])
            mat[idx] = dij
        end
    end
    dropzeros!(mat)
    return param
end

function dufal_mix(bondvol::AssocParam,::Nothing)
    length(bondvol.values.values) == 0 && return deepcopy(bondvol)

    param = assoc_extend(bondvol)
    mat = param.values
    for (idx,(i,j),(a,b)) in indices(mat)
        i == j && continue
        dij = mix_mean3(mat[i,i][a,b],mat[j,j][a,b])
        if !iszero(dij) && iszero(mat[idx])
            mat[idx] = dij
        end
    end
    dropzeros!(mat)
    return param
end


function epsilon_assoc_mix(epsilon_assoc::AssocParam)
    length(epsilon_assoc.values.values) == 0 && return deepcopy(epsilon_assoc)

    param = assoc_extend(epsilon_assoc)
    mat = param.values
    for (idx,(i,j),(a,b)) in indices(mat)
        i == j && continue
        dij = (mat[i,i][a,b] + mat[j,j][a,b])/2
        if !iszero(dij) && iszero(mat[idx])
            mat[idx] = dij
        end
    end
    dropzeros!(mat)
    return param
end

function bondvol_mix(bondvol::AssocParam,σ)
    length(bondvol.values.values) == 0 && return deepcopy(bondvol)

    param = assoc_extend(bondvol)
    mat = param.values
    for (idx,(i,j),(a,b)) in indices(mat)
        i == j && continue
        dij = sqrt(mat[i,i][a,b]*mat[j,j][a,b])*(sqrt(σ[i,i]*σ[j,j])/σ[i,j])^3
        if !iszero(dij) && iszero(mat[idx])
            mat[idx] = dij
        end
    end
    dropzeros!(mat)
    return param
end

function zero_mix(assocparam::AssocParam)
    length(assocparam.values.values) == 0 && return deepcopy(assocparam)
    param = assoc_extend(assocparam)
    mat = param.values

    #fill with sentinel values
    SENTINEL = -124
    for (idx,(i,j),(a,b)) in indices(mat)
        i == j && continue
        di = mat[i,i][a,b]
        dj = mat[j,j][a,b]
        dij = max(di*dj,zero(di*dj))
        if !iszero(dij) && iszero(mat[idx])
            mat.values[idx] = SENTINEL
        end
    end
    dropzeros!(mat)
    #refill with zeros.
    for i in eachindex(mat.values)
        if mat[i] == SENTINEL
            mat[i] = 0
        end
    end

    return param
end

function assoc_mix(bondvol,epsilon_assoc,sigma,assoc_options::AssocOptions)
    combining = assoc_options.combining
    if combining == :nocombining
        return bondvol,epsilon_assoc
    elseif combining in (:elliott_runtime,:esd_runtime)
        #return bondvol,epsilon_assoc
        return zero_mix(bondvol),zero_mix(epsilon_assoc)
    elseif combining in (:elliott,:esd)
        return bondvol_mix(bondvol,sigma),epsilon_assoc_mix(epsilon_assoc)
    elseif combining == :cr1
        return bondvol_mix(bondvol,nothing),epsilon_assoc_mix(epsilon_assoc)
    elseif combining in (:dufal,:mie15)
        return dufal_mix(bondvol,nothing),epsilon_assoc_mix(epsilon_assoc)
    else
        throw(error("incorrect combining argument ",error_color(string(combining))," passed to AssocOptions."))
    end
end

function assoc_mix!(data,components)
    assoc_options = data["assoc_options"]
    if haskey(data,"bondvol") && haskey(data,"epsilon_assoc")
        bondvol = data["bondvol"]
        epsilon_assoc = data["epsilon_assoc"]
        sigma = get(data,"sigma",nothing)
        bondvol, epsilon_assoc = assoc_mix(bondvol,epsilon_assoc,sigma,assoc_options)
        data["bondvol"] = bondvol
        data["epsilon_assoc"] = epsilon_assoc
    else
        x1 = get!(data,"bondvol") do
            AssocParam("bondvol",components)
        end
        x2 = get!(data,"epsilon_assoc") do
            AssocParam("epsilon_assoc",components)
        end
    end
    return data
end
