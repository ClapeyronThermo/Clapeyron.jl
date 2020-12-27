abstract type EoS end
abstract type SAFT <: EoS end
abstract type Cubic <: EoS end
abstract type Ideal <: EoS end

abstract type PCSAFTFamily <: SAFT end
abstract type CKSAFTFamily <: SAFT end
abstract type sPCSAFTFamily <: SAFT end
abstract type ogSAFTFamily <: SAFT end
abstract type SAFTVRSWFamily <: SAFT end
abstract type SAFTVRMieFamily <: SAFT end
abstract type SAFTVRQMieFamily <: SAFT end
abstract type SAFTgammaMieFamily <: SAFT end
abstract type softSAFTFamily <: SAFT end

abstract type vdWFamily <: Cubic end
abstract type RKFamily <: Cubic end
abstract type SRKFamily <: Cubic end
abstract type PRFamily <: Cubic end
abstract type CPAFamily <: Cubic end

#to allow broadcasting without problems
Base.broadcastable(model::EoS) = Ref(model)


struct SAFTVRMie <: SAFTVRMieFamily
    components::Array{Set{String},1}
    params::SAFTVRMieParams
    ideal::Ideal
end

struct SAFTVRQMie <: SAFTVRQMieFamily
    components::Array{Set{String},1}
    params::SAFTVRQMieParams
    ideal::Ideal
end

struct SAFTVRSW <: SAFTVRSWFamily
    components::Array{Set{String},1}
    params::SAFTVRSWParams
    ideal::Ideal
end

struct CKSAFT <: CKSAFTFamily
    components::Array{Set{String},1}
    params::CKSAFTParams
    ideal::Ideal
end

struct PCSAFT <: PCSAFTFamily
    components::Array{Set{String},1}
    sites::Array{String,1}
    params::PCSAFTParams
    ideal::Ideal
end

struct sPCSAFT <: sPCSAFTFamily
    components::Array{Set{String},1}
    params::sPCSAFTParams
    ideal::Ideal
end

struct ogSAFT <: ogSAFTFamily
    components::Array{Set{String},1}
    params::ogSAFTParams
    ideal::Ideal
end

struct softSAFT <: softSAFTFamily
    components::Array{Set{String},1}
    params::softSAFTParams
    ideal::Ideal
end

struct SAFTgammaMie <: SAFTgammaMieFamily
    components::Array{Set{String},1}
    groups::Array{Set{String},1}
    group_multiplicities::Dict{Set{String},DefaultDict{Set{String},Int64,Int64}}
    params::SAFTgammaMieParams
    ideal::Ideal
end

struct vdW <: vdWFamily
    components::Array{Set{String},1}
    params::vdWParams
    ideal::Ideal
end

struct RK <: RKFamily
    components::Array{Set{String},1}
    params::RKParams
    ideal::Ideal
end

struct SRK <: SRKFamily
    components::Array{Set{String},1}
    params::SRKParams
    ideal::Ideal
end

struct PR <: PRFamily
    components::Array{Set{String},1}
    params::PRParams
    ideal::Ideal
end

struct CPA <: CPAFamily
    components::Array{Set{String},1}
    params::CPAParams
    ideal::Ideal
end

# Ideal #

struct Monomer <: Ideal
    components::Array{Set{String},1}
    params::MonomerParams
end

struct Wilhoit <: Ideal
    params::WilhoitParams
end

struct NASA <: Ideal
    params::NASAParams
end

struct Walker <: Ideal
    components::Array{Set{String},1}
    params::WalkerParams
end

struct Reid <: Ideal
    components::Array{Set{String},1}
    params::ReidParams
end

struct Basic <: Ideal
    components::Array{Set{String},1}
end
