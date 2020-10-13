abstract type EoS end
abstract type SAFT <: EoS end

abstract type PCSAFTFamily <: SAFT end
abstract type sPCSAFTFamily <: SAFT end
abstract type ogSAFTFamily <: SAFT end
abstract type SAFTVRMieFamily <: SAFT end
abstract type SAFTgammaMieFamily <: SAFT end

struct SAFTVRMie <: SAFTVRMieFamily
    components::NamedArray{Float64,1}
    params::SAFTVRMieParams
end

struct PCSAFT <: PCSAFTFamily
    components::NamedArray{Float64,1}
    params::PCSAFTParams
end

struct sPCSAFT <: sPCSAFTFamily
    components::NamedArray{Float64,1}
    params::sPCSAFTParams
end

struct ogSAFT <: ogSAFTFamily
    components::NamedArray{Float64,1}
    params::ogSAFTParams
end

struct SAFTgammaMie <: SAFTgammaMieFamily
    components::Array{AbstractString,1}
    groups::Dict{AbstractString,Int64}
    params::SAFTgammaMieParams
end
