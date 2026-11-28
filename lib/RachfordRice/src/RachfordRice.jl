module RachfordRice

include("misc.jl")
include("struct.jl")
include("rr.jl")
include("check_result.jl")

"""
    vapour_compositions(K,z,β)
    vapour_compositions!(y,K,z,β)
    vapour_compositions(RR::RachfordRiceResidual,β)
    vapour_compositions!(y,RR::RachfordRiceResidual,β)

Returns the gas phase composition, given K-values `K`, the initial molar composition `z` and the molar vapor fraction ´β´.

Each gas phase composition is calculated acording to:

xvᵢ = kᵢzᵢ/(1+ β(kᵢ - 1))

`vapour_compositions!` does the same, but it overwrites the input `y`.
"""
function vapour_compositions(K, z, βv, βl=1 - βv)
    TT = Base.promote_eltype(K, z, βv, βl)
    vapour_compositions!(similar(z, TT), K, z, βv, βl)
end

function vapour_compositions(RR::RachfordRiceResidual, β_result::RachfordRiceResult)
    βv, βl = β_result.V, β_result.L
    vapour_compositions(RR, RR.z, βv, βl)
end

function vapour_compositions!(y, K, z, βv, βl=1 - βv)
    _1 = one(eltype(y))
    n = sum(z)
    for i in eachindex(y)
        Ki, zi = K[i], z[i]
        wi = zi/n
        if isinf(Ki)
            y[i] = wi/(βv)
        else
            detKi = inv_fma(Ki, βv, βl)
            y[i] = detKi*wi*Ki
        end
    end
    return y
end

vapour_compositions(RR::RachfordRiceResidual, βv, βl=1 - βv) = vapour_compositions(RR, RR.z, βv, βl)
vapour_compositions!(y, RR::RachfordRiceResidual, βv, βl=1 - βv) = vapour_compositions!(y, RR, RR.z, βv, βl)

"""
    liquid_compositions(K,z,β)
    liquid_compositions!(x,K,z,β)
    liquid_compositions(RR::RachfordRiceResidual,β::Union{Number,RachfordRiceResult})
    liquid_compositions!(x,RR::RachfordRiceResidual,β::Union{Number,RachfordRiceResult})

Returns the liquid phase composition, given k-values `K`, the initial molar composition `z` and the molar vapor fraction ´β´.

Each gas phase composition is calculated acording to:

    xlᵢ = zᵢ/(1+ β(kᵢ-1))

`liquid_compositions!` does the same, but it overwrites the input `x`.
"""
function liquid_compositions(K, z, βv, βl=1 - βv)
    TT = Base.promote_eltype(K, z, βv, βl)
    liquid_compositions!(similar(z, TT), K, z, βv, βl)
end

function liquid_compositions(RR::RachfordRiceResidual, β_result::RachfordRiceResult)
    βv, βl = β_result.V, β_result.L
    liquid_compositions(RR, RR.z, βv, βl)
end

function liquid_compositions!(x, K, z, βv, βl=1 - βv)
    _1 = one(eltype(x))
    _0 = zero(_1)
    n = sum(z)
    for i in eachindex(x)
        Ki, zi = K[i], z[i]
        if isinf(Ki)
            x[i] = _0
        else
            wi = zi/n
            detKi = inv_fma(Ki, βv, βl)
            x[i] = detKi*wi
        end
    end
    return x
end

liquid_compositions(RR::RachfordRiceResidual, βv, βl=1 - βv) = liquid_compositions(RR, RR.z, βv, βl)
liquid_compositions!(y, RR::RachfordRiceResidual, βv, βl=1 - βv) = liquid_compositions!(y, RR, RR.z, βv, βl)

"""
    liquid_compositions(K,z,β)
    liquid_compositions!(x,K,z,β)

Returns the liquid phase composition, given k-values `K`, the initial molar composition `z` and the molar vapor fraction ´β´.

Each gas phase composition is calculated acording to:

    xlᵢ = zᵢ/(1+ β(kᵢ-1))

`liquid_compositions!` does the same, but it overwrites the input `x`.
"""
liquid_compositions!

"""
    vapour_compositions(K,z,β)
    vapour_compositions!(y,K,z,β)
    vapour_compositions(RR::RachfordRiceResidual,β::Union{Number,RachfordRiceResult})
    vapour_compositions!(y,RR::RachfordRiceResidual,β::Union{Number,RachfordRiceResult})

Returns the gas phase composition, given K-values `K`, the initial molar composition `z` and the molar vapor fraction ´β´.

Each gas phase composition is calculated acording to:

xvᵢ = kᵢzᵢ/(1+ β(kᵢ - 1))

`vapour_compositions!` does the same, but it overwrites the input `y`.
"""
vapour_compositions!

export rachfordrice
export liquid_compositions, vapour_compositions, liquid_compositions!, vapour_compositions!
export RachfordRiceResidual, RachfordRiceResult
export RachfordRiceStatus, RRLiquid, RRVapour, RREq, RRFailure, RRTrivial

@public residual
@public residual_first_derivative
@public residual_first_second_derivative
@public residual_first_second_derivative_error

@public rr_limits
@public rr_feasible
@public rr_status

end #module
