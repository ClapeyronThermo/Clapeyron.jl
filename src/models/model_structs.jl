abstract type EoS end
abstract type SAFT <: EoS end

abstract type PCSAFTFamily <: SAFT end
abstract type sPCSAFTFamily <: SAFT end
abstract type ogSAFTFamily <: SAFT end
abstract type SAFTVRMieFamily <: SAFT end

struct SAFTVRMie <: SAFTVRMieFamily; 
    components
    params::SAFTVRMieParams
end

function SAFTVRMie(components::T1) where T1 <: AbstractVector{T2} where T2 <: AbstractString
    set_components = [Set([components[i]]) for i in 1:length(components)]
    raw_params = retrieveparams(components, "SAFTVRMie"; kwargs...)
    model = SAFTVRMie(set_components, create_SAFTVRMieParams(raw_params))
end

struct PCSAFT <: PCSAFTFamily
    components
    params::PCSAFTParams
end

function PCSAFT(components::T1) where T1 <: AbstractVector{T2} where T2 <: AbstractString
    set_components = [Set([components[i]]) for i in 1:length(components)]
    raw_params = retrieveparams(components, "PCSAFT"; kwargs...)
    model = PCSAFT(set_components, create_PCSAFTParams(raw_params))
end

struct sPCSAFT <: sPCSAFTFamily
    components
    params::sPCSAFTParams
end

function sPCSAFT(components::T1) where T1 <: AbstractVector{T2} where T2 <: AbstractString
    set_components = [Set([components[i]]) for i in 1:length(components)]
    raw_params = retrieveparams(components, "sPCSAFT"; kwargs...)
    model = sPCSAFT(set_components, create_sPCSAFTParams(raw_params))
end

struct ogSAFT <: ogSAFTFamily
    components
    params::ogSAFTParams
end

function ogSAFT(components::T1) where T1 <: AbstractVector{T2} where T2 <: AbstractString
    set_components = [Set([components[i]]) for i in 1:length(components)]
    raw_params = retrieveparams(components, "ogSAFT"; kwargs...)
    model = ogSAFT(set_components, create_ogSAFTParams(raw_params))
end
