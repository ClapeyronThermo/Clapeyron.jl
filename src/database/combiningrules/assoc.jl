function assoc_extend(param::AssocParam)
    length(param.values.values) == 0 && return param
    _4dmatrix = assoc_extend(param.values,param.sites)
    return AssocParam(param.name,param.components,_4dmatrix,param.sites,param.sourcecsvs,param.sources)
end

function assoc_extend(mat::Compressed4DMatrix,sites)
    length(mat.values) == 0 && return mat
    #c = length(param)
    comps = length(sites)
    vals,c,s = mat.values,mat.outer_indices,mat.inner_indices
    idx = Vector{NTuple{4,Int}}(undef,0)
    for i in 1:comps
        for j in 1:i #include diagonal
            la = length(sites[i])
            lb = length(sites[j])
            if (la != 0) && (lb !=0) #delete empty interactions
                for a in 1:la
                    #when i == j, we are in association site pairs of a single component. those are symmetrical.
                    #the same cannot be said of intercomponent association pairs. that is x[i,j][a,b] could be different that x[j,i][a,b]
                    start = ifelse(i == j,a,1)
                    for b in start:lb
                        push!(idx,(i,j,a,b))
                    end
                end
            end
        end
    end
    sort!(idx)
    extended_vals = zeros(eltype(vals),length(idx))
    for (k,(i,j,a,b)) in enumerate(idx)
        extended_vals[k] = mat[i,j][a,b]
    end
    extended_outer_indices = [(c[1],c[2]) for c ∈ idx]
    extended_inner_indices = [(c[3],c[4]) for c ∈ idx]
    #unsafe constructor
    return Compressed4DMatrix(extended_vals,extended_outer_indices,extended_inner_indices,true)
end

bondvol_mix(bondvol::AssocParam) = bondvol_mix(bondvol,nothing)

function bondvol_mix(bondvol::AssocParam,::Nothing)
    length(bondvol.values.values) == 0 && return deepcopy(bondvol)
    param = assoc_extend(bondvol)
    mat = param.values
    for (idx,(i,j),(a,b)) in indices(mat)
        if iszero(mat.values[idx])
            mat.values[idx] = sqrt(mat[i,i][a,b]*mat[j,j][a,b])
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
        if iszero(mat.values[idx])
            mat.values[idx] = (mat[i,i][a,b] + mat[j,j][a,b])/2
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
        if iszero(mat.values[idx])
            mat.values[idx] = sqrt(mat[i,i][a,b]*mat[j,j][a,b])*(sqrt(σ[i,i]*σ[j,j])/σ[i,j])^3
        end
    end
    dropzeros!(mat)
    return param
end

function assoc_mix(bondvol,epsilon_assoc,sigma,assoc_options::AssocOptions)
    combining = assoc_options.combining
    if combining == :nocombining
        return bondvol,epsilon_assoc
    elseif combining in (:elliott_runtime,:esd_runtime)
        if assoc_options.dense
            return bondvol,epsilon_assoc
        else
            throw(error("cannot use sparse solver with :elliot_runtime combining rule"))
        end
    elseif combining in (:elliott,:esd)
        return bondvol_mix(bondvol,sigma),epsilon_assoc_mix(epsilon_assoc)
    elseif combining == :cr1
        return bondvol_mix(bondvol),epsilon_assoc_mix(epsilon_assoc)
    else
        throw(error("incorrect combining argument ",error_color(string(combining))," passed to AssocOptions."))
    end
end

function assoc_mix!(data,components,assoc_options = AssocOptions())
    if haskey(data,"bondvol") && haskey(data,"epsilon_assoc")  
        bondvol = data["bondvol"]
        epsilon_assoc = data["epsilon_assoc"]
        sigma = get(data,"sigma",nothing)
        bondvol, epsilon_assoc = assoc_mix(bondvol,epsilon_assoc,sigma,assoc_options)
        data["bondvol"] = bondvol
        data["epsilon_assoc"] = epsilon_assoc
    else
        data["bondvol"] = AssocParam("bondvol",components)
        data["epsilon_assoc"] = AssocParam("epsilon_assoc",components)
    end
    return data
end