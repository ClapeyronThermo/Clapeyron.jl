include("../../../../lib/RachfordRice/src/RachfordRice.jl")

using Clapeyron.RachfordRice
using Clapeyron.RachfordRice: rr_feasible

function material_balance_rr_converged(w,z,β::Number,n = sum(z),ztol = sqrt(eps(eltype(β))))
    βx = SVector(1 - β,β)
    return material_balance_rr_converged(w,z,βx,n,ztol)
end

function material_balance_rr_converged(w,z,β::AbstractVector,n = sum(z),ztol = sqrt(eps(eltype(β))))
    rk = zero(Base.promote_eltype(w[1],z,β))
    np = length(β)
    #material balance test from https://github.com/WhitsonAS/Rachford-Rice-Contest
    for i in eachindex(z)
        rk1 = zero(rk)
        rk2 = zero(rk)
        zi = z[i]/n
        for j in 1:np
            Vi = β[j]*w[j][i]
            rk1 +=  Vi
            rk2 += abs(Vi)
        end
        if !iszero(zi)
            rk = max(rk,abs(rk1 - zi)/abs(rk2 + zi))
        end
    end
    return rk <= ztol
end


"""
    rr_feasibility_projection!(K, z, bubble::Bool, verbose = false) -> β_edge

rescale all `K` values by `c`, with the scaled K vector presenting bubble or dew equilibria:

    - bubble: c = Σzᵢ / ΣzᵢKᵢ (finite)
    - dew: c = Σzᵢ/Kzᵢ (finite)

Returns the vapour fraction of the rescaled system, (0 or 1 for normal systems, 0 < β < 1 for system that have non-volatile or non-condensable components)
"""
function rr_feasibility_projection!(K, z, bubble; verbose = false)
    if bubble
        K0,β = bubble_scaling_restriction(K,z)
        verbose && @info "maximum(K) < 1: projecting K-values to bubble point"
        if 0 < β < 1
            verbose && @info "supposing β = $β (bubble initialization with non-condensable components)" 
        else
            verbose && @info "supposing β = 0 (bubble initialization)"
        end
    else
        K0,β = dew_scaling_restriction(K,z)
        if 0 < β < 1
            verbose && @info "supposing β = $β (dew initialization with non-volatile components)" 
        else
            verbose && @info "supposing β = 0 (dew initialization)"
        end
    end
    for i in eachindex(K)
        K[i] = K[i]*K0
    end
    verbose && @info "scaled K: $K"
    return β
end

function bubble_scaling_restriction(K, z)
    T = Base.promote_eltype(K, z)
    ∑z = zero(T); ∑K = zero(T); Z∞ = zero(T)
    @inbounds for i in eachindex(K)
        Ki, zi = K[i], z[i]
        iszero(zi) && continue
        if isinf(Ki)
            Z∞ += zi
        else
            ∑z += zi
            ∑K += zi * Ki
        end
    end
    ∑zt = ∑z + Z∞
    return (∑z / ∑K, Z∞/∑zt)
end

function dew_scaling_restriction(K, z)
    T = Base.promote_eltype(K, z)
    ∑invK = zero(T); ∑z = zero(T); Z₀ = zero(T)
    @inbounds for i in eachindex(K)
        Ki, zi = K[i], z[i]
        iszero(zi) && continue
        if iszero(Ki)
            Z₀ += zi
        else
            ∑z += zi
            ∑invK += zi / Ki
        end
    end
    ∑zt = ∑z + Z₀
    return (∑invK / ∑z, one(T) - Z₀/∑zt)
end