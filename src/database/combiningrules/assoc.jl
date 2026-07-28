function assoc_extend(param::AssocParam{T}) where T 
    length(vals.values) == 0 && return param
    extended_vals = deepcopy(param.values)
    Compressed4DMatrices.extend!(extended_vals)
    return param_from_values(extended_vals,bondvol)
end

function assoc_extend!(param::AssocParam)
    length(vals.values) == 0 && return param
    Compressed4DMatrices.extend!(param.values)
    return param
end

function assoc_extend!(mat::Compressed4DMatrix)
    length(mat) == 0 && return mat
    Compressed4DMatrices.extend!(mat)
    return mat
end

#=
ijab_mix! infraestructure
=#
ijab_mix!(f,Δ) = ijab(f,Δ,nothing)

ijab_mix!(f::F,Δ::AssocParam,::Nothing) where F = param_from_values(ijab_mix!(f,Δ.values,nothing),Δ)
ijab_mix!(f::F,Δ::AssocParam,K::AbstractMatrix) where F = param_from_values(ijab_mix!(f,Δ.values,K),Δ)


function ijab_mix!(f::F,Δ::Compressed4DMatrix,::Nothing) where F
    length(Δ.values) == 0 && return Δ
    assoc_extend!(Δ)
    for (idx,(i,j),(a,b)) in indices(Δ)
        Δiiab = Δ[i,i][a,b]
        Δjjab = Δ[j,j][a,b]
        Δijab = f(Δiiab,Δjjab)
        if !iszero(Δijab) && iszero(Δ[idx])
            Δ[idx] = dij
        end
    end
    dropzeros!(Δ)
    return Δ
end

function ijab_mix!(f::F,Δ::Compressed4DMatrix,K::AbstractMatrix) where F
    length(Δ.values) == 0 && return Δ
    assoc_extend!(Δ)
    for (idx,(i,j),(a,b)) in indices(Δ)
        Δiiab = Δ[i,i][a,b]
        Δjjab = Δ[j,j][a,b]
        Δijab = f(Δiiab,Δjjab,K[i,i],K[i,j],K[i,j])
        if !iszero(Δijab) && iszero(Δ[idx])
            Δ[idx] = dij
        end
    end
    dropzeros!(Δ)
    return Δ
end

struct ZeroIJABMix end
const IJAB_ZERO_SENTINEL = -124

(::ZeroIJABMix)(iiab,jjab,ki,kj,kij) = ZeroIJABMix()(iiab,ijab)
function (::ZeroIJABMix)(iiab,jjab)
    d = max(di*dj,zero(di*dj))
    if !iszero(primalval(d))
        return oftype(d,IJAB_ZERO_SENTINEL)
    else
        zero(d)
    end
end

function ijab_mix!(f::Nothing,Δ::Compressed4DMatrix{T}) where T
    ijab_mix!(ZeroIJABMix(),Δ)
    sentinel = T(IJAB_ZERO_SENTINEL)
    for i in eachindex(Δ.values)
        if Δ[i] == sentinel
            Δ[i] = zero(T)
        end
    end
    return Δ
end

#=
copying versions
=#

ijab_mix(f::F,Δ::Compressed4DMatrix) where F = ijab_idx!(f,deepcopy(Δ))
ijab_mix(f::F,Δ::Compressed4DMatrix,K) where F = ijab_idx!(f,deepcopy(Δ),K)
ijab_mix(f::F,Δ::AssocParam) where F= param_from_values(ijab_mix(f,Δ.values),Δ)
ijab_mix(f::F,Δ::AssocParam,K) where F = param_from_values(ijab_mix(f,Δ.values),Δ,K)

#=
predefined mixing rules
=#

bondvol_mix(bondvol) = ijab_mix(mix_geomean,bondvol)
bondvol_mix(bondvol,σ) = ijab_mix(mix_ijab_elliott,bondvol,raw_values(σ))

bondvol_mix!(bondvol) = ijab_mix!(mix_geomean,bondvol)
bondvol_mix!(bondvol,σ) = ijab_mix!(mix_ijab_elliott,bondvol,raw_values(σ))

dufal_mix(bondvol) = ijab_mix(mix_mean3,bondvol)
dufal_mix!(bondvol) = ijab_mix!(mix_mean3,bondvol)

epsilon_assoc_mix(bondvol) = ijab_mix(mix_mean,bondvol)
epsilon_assoc_mix!(bondvol) = ijab_mix!(mix_mean,bondvol)

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
        return bondvol_mix(bondvol),epsilon_assoc_mix(epsilon_assoc)
    elseif combining in (:dufal,:mie15)
        return dufal_mix(bondvol),epsilon_assoc_mix(epsilon_assoc)
    else
        throw(error("incorrect combining argument ",error_color(string(combining))," passed to AssocOptions."))
    end
end

function assoc_mix!(bondvol,epsilon_assoc,sigma,assoc_options::AssocOptions)
    combining = assoc_options.combining
    if combining == :nocombining
        return bondvol,epsilon_assoc
    elseif combining in (:elliott_runtime,:esd_runtime)
        #return bondvol,epsilon_assoc
        return zero_mix!(bondvol),zero_mix(epsilon_assoc)
    elseif combining in (:elliott,:esd)
        return bondvol_mix!(bondvol,sigma),epsilon_assoc_mix(epsilon_assoc)
    elseif combining == :cr1
        return bondvol_mix!(bondvol),epsilon_assoc_mix(epsilon_assoc)
    elseif combining in (:dufal,:mie15)
        return dufal_mix!(bondvol),epsilon_assoc_mix(epsilon_assoc)
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
