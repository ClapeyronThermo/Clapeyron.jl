function assoc_extend(param::AssocParam{T}) where T 
    length(param.values) == 0 && return param
    extended_vals = deepcopy(param.values)
    Compressed4DMatrices.extend!(extended_vals)
    return param_from_values(extended_vals,param)
end

function assoc_extend!(param::AssocParam)
    length(param.values) == 0 && return param
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
"""
    ijab_mix!(f,Δ::AssocParam)
    ijab_mix!(f,Δ::AssocParam,K::AbstractMatrix)
    ijab_mix!(f,Δ::Compressed4DMatrix,K = nothing)
    ijab_mix!(nothing,Δ)

Inplace combining rule for association (site-site) parameters, such as bond volumes (`κᵢⱼₐᵦ`) or association energies (`εᵢⱼₐᵦ`).

Given the pure-component (diagonal, `i == i`) association values already stored in `Δ`, `ijab_mix!` fills in the *unset* cross-association entries (`i != j`) using a two-body combining rule `f`, and overwrites `Δ` with the result. Only entries that are currently zero (unset) are populated; explicitly-provided cross parameters are left untouched.

For each stored site pair `(a,b)` shared between components `i` and `j`, the new cross term is calculated as:
```
Δᵢⱼₐᵦ = f(Δᵢᵢₐᵦ, Δⱼⱼₐᵦ)          # if K === nothing
Δᵢⱼₐᵦ = f(Δᵢᵢₐᵦ, Δⱼⱼₐᵦ, Kᵢᵢ, Kᵢⱼ, Kᵢⱼ) # if K is a matrix
```
Where `f` is a two-argument (or five-argument, when a `K` matrix is supplied) 'combining' function, analogous in spirit to the ones used by [`kij_mix`](@ref)/[`pair_mix`](@ref) but operating on association-site data instead of plain pair data.

If `f == nothing`, `ijab_mix!` will only add the corresponding mixing entries with the value of zero.

Mutates and returns `Δ`. If `Δ` (or its underlying `Compressed4DMatrix`) has no entries, it is returned unchanged.

See also the non-mutating counterpart [`ijab_mix`](@ref).
"""
ijab_mix!(f,Δ) = ijab_mix!(f,Δ,nothing)

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
            Δ[idx] = Δijab
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
        Δijab = f(Δiiab,Δjjab,K[i,i],K[j,j],K[i,j])
        if !iszero(Δijab) && iszero(Δ[idx])
            Δ[idx] = Δijab
        end
    end
    dropzeros!(Δ)
    return Δ
end

struct ZeroIJABMix end
const IJAB_ZERO_SENTINEL = -124

(::ZeroIJABMix)(iiab,jjab,ki,kj,kij) = ZeroIJABMix()(iiab,jjab)
function (::ZeroIJABMix)(iiab,jjab)
    Δ = max(iiab*jjab,zero(iiab*jjabj))
    if !iszero(primalval(Δ))
        return oftype(Δ,IJAB_ZERO_SENTINEL)
    else
        zero(Δ)
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

"""
    ijab_mix(f,Δ::AssocParam)
    ijab_mix(f,Δ::AssocParam,K)
    ijab_mix(f,Δ::Compressed4DMatrix,K = nothing)

Non-mutating version of [`ijab_mix!`](@ref).

Combines the cross-association entries of `Δ` using the two-body (or, with `K`, five-argument) rule `f`, following exactly the same rules and formulas as [`ijab_mix!`](@ref), but operates on a copy of `Δ` (or of its underlying `Compressed4DMatrix`) so that the original parameter is left untouched.

Returns a new `AssocParam`/`Compressed4DMatrix` of the same shape as `Δ`, with previously-unset cross terms populated by `f` and existing entries left as they were.
"""
ijab_mix(f::F,Δ::Compressed4DMatrix) where F = ijab_mix!(f,deepcopy(Δ))
ijab_mix(f::F,Δ::Compressed4DMatrix,K) where F = ijab_mix!(f,deepcopy(Δ),K)
ijab_mix(f::F,Δ::AssocParam) where F= param_from_values(ijab_mix(f,Δ.values),Δ)
ijab_mix(f::F,Δ::AssocParam,K) where F = param_from_values(ijab_mix(f,Δ.values,K),Δ)

#=
predefined mixing rules
=#

"""
    bondvol_mix(bondvol::AssocParam)
    bondvol_mix(bondvol::AssocParam,σ)

Combining rule for cross-association bond volumes (`κᵢⱼₐᵦ`), used by the `:cr1` and `:elliott`/`:esd` combining rules.

- `bondvol_mix(bondvol)` (CR-1 combining rule) fills unset cross terms with the geometric mean of the pure-component bond volumes:
```
κᵢⱼₐᵦ = √(κᵢᵢₐᵦ * κⱼⱼₐᵦ)
```

- `bondvol_mix(bondvol,σ)` (Elliott/ESD-style combining rule) additionally uses the segment/site diameters `σ` to correct the geometric mean for size differences between components:
```
κᵢⱼₐᵦ = √(κᵢᵢₐᵦ * κⱼⱼₐᵦ * σᵢ^3 * σⱼ^3) / σᵢⱼ^3
```

See [`bondvol_mix!`](@ref) for the inplace version.
"""
bondvol_mix(bondvol) = ijab_mix(mix_geomean,bondvol)
bondvol_mix(bondvol,σ) = ijab_mix(mix_ijab_elliott,bondvol,raw_values(σ))

"""
    bondvol_mix!(bondvol::AssocParam)
    bondvol_mix!(bondvol::AssocParam,σ)

Inplace version of [`bondvol_mix`](@ref). Mutates `bondvol`, filling in unset cross-association bond volumes in place using either the plain geometric mean, one-argument form) or the size-corrected Elliott/ESD rule, two-argument form using segment diameters `σ`). 
See [`bondvol_mix`](@ref) for the exact formulas.
"""
bondvol_mix!(bondvol) = ijab_mix!(mix_geomean,bondvol)
bondvol_mix!(bondvol,σ) = ijab_mix!(mix_ijab_elliott,bondvol,raw_values(σ))

"""
    dufal_mix(bondvol::AssocParam)

Combining rule for cross-association bond volumes (`κᵢⱼₐᵦ`) following the Dufal (SAFT-VR Mie / Mie 15) mixing rule, used by the `:dufal`/`:mie15` combining rules.

Fills unset cross terms with the cubic-mean-radius rule:
```
κᵢⱼₐᵦ = (0.5*(∛κᵢᵢₐᵦ + ∛κⱼⱼₐᵦ))^3
```

See [`dufal_mix!`](@ref) for the inplace version.
"""
dufal_mix(bondvol) = ijab_mix(mix_mean3,bondvol)

"""
    dufal_mix!(bondvol::AssocParam)

Inplace version of [`dufal_mix`](@ref). 
Mutates `bondvol`, filling in unset cross-association bond volumes in place using the Dufal cubic-mean-radius rule. 
See [`dufal_mix`](@ref) for the exact formula.
"""
dufal_mix!(bondvol) = ijab_mix!(mix_mean3,bondvol)

"""
    epsilon_assoc_mix(epsilon_assoc::AssocParam)

Combining rule for cross-association energies (`εᵢⱼₐᵦ`), used across all combining rules that mix association energy (`:cr1`, `:elliott`/`:esd`, `:dufal`/`:mie15`.

Fills unset cross terms with the arithmetic mean of the pure-component association energies:
```
εᵢⱼₐᵦ = 0.5*(εᵢᵢₐᵦ + εⱼⱼₐᵦ)
```

See [`epsilon_assoc_mix!`](@ref) for the inplace version.

"""
epsilon_assoc_mix(epsilon_assoc) = ijab_mix(mix_mean,epsilon_assoc)

"""
    epsilon_assoc_mix!(epsilon_assoc::AssocParam)

Inplace version of [`epsilon_assoc_mix`](@ref). 
Mutates its argument, filling in unset cross-association energies in place using the arithmetic mean rule. 
See [`epsilon_assoc_mix`](@ref) for the exact formula.
"""
epsilon_assoc_mix!(bondvol) = ijab_mix!(mix_mean,epsilon_assoc)

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

assoc_mix(bondvol,epsilon_assoc,assoc_options::AssocOptions) = assoc_mix(bondvol,epsilon_assoc,nothing,assoc_options)
assoc_mix!(bondvol,epsilon_assoc,assoc_options::AssocOptions) = assoc_mix!(bondvol,epsilon_assoc,nothing,assoc_options)

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