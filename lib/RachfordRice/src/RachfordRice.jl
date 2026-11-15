module RachfordRice

include("misc.jl")
include("struct.jl")
include("rr.jl")

"""
    vapour_compositions(K,z,β)
    vapour_compositions!(y,K,z,β)
    vapour_compositions(RR::RachfordRiceResidual,β)
    vapour_compositions!(y,RR::RachfordRiceResidual,β)

Returns the gas phase composition, given K-values `K`, the initial molar composition `z` and the molar vapor fraction ´β´.

Each gas phase composition is calculated acording to:

xvᵢ = kᵢzᵢ/(1+ β(kᵢ - 1))
"""
function vapour_compositions(K,z,β)
    TT = Base.promote_eltype(K,z,β)
    vapour_compositions!(similar(z,TT),K,z,β)
end

function vapour_compositions!(y,K,z,β)
    _1 = one(eltype(y))
    _β = _1*β
    for i in eachindex(y)
        Ki,zi = K[i],z[i]
        y[i] = ifelse(isinf(Ki), zi/_β, Ki*zi/(_1 + _β*(Ki -_1)))
    end
    return y
end

vapour_compositions(RR::RachfordRiceResidual,β) = vapour_compositions(RR,RR.z,β)
vapour_compositions!(y,RR::RachfordRiceResidual,β) = vapour_compositions!(y,RR,RR.z,β)

"""
    liquid_compositions(K,z,β)
    liquid_compositions!(x,K,z,β)

Returns the liquid phase composition, given k-values `K`, the initial molar composition `z` and the molar vapor fraction ´β´.

Each gas phase composition is calculated acording to:

    xlᵢ = zᵢ/(1+ β(kᵢ-1))
"""
function liquid_compositions(K,z,β)
    TT = Base.promote_eltype(K,z,β)
    liquid_compositions!(similar(z,TT),K,z,β)
end

function liquid_compositions!(x,K,z,β)
    _1 = one(eltype(x))
    _β = _1*β
    for i in eachindex(x)
        Ki,zi = K[i],z[i]
        x[i] = ifelse(isinf(Ki), zero(_1), zi/(_1 + _β*(Ki -_1)))
    end
    return x
end

liquid_compositions(RR::RachfordRiceResidual,β) = liquid_compositions(RR,RR.z,β)
liquid_compositions!(y,RR::RachfordRiceResidual,β) = liquid_compositions!(y,RR,RR.z,β)

export rachfordrice
export liquid_compositions,vapour_compositions,liquid_compositions!,vapour_compositions!
export RachfordRiceResidual, RachfordRiceResult, vapour_fraction, liquid_fraction
export RachfordRiceStatus, RRLiquid, RRVapour, RREq, RRFailure, RRTrivial

@public residual
@public residual_first_derivative
@public residual_first_second_derivative
@public residual_first_second_derivative_error

end #module
