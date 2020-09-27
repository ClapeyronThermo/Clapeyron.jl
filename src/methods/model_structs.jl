abstract type EoS end
abstract type SAFT <: EoS end

abstract type PCSAFTFamily <: SAFT end
abstract type sPCSAFTFamily <: SAFT end
abstract type ogSAFTFamily <: SAFT end
abstract type SAFTVRMieFamily <: SAFT end

struct SAFTVRMie <: SAFTVRMieFamily; components; parameters::SAFTVRMieParams end
struct PCSAFT <: PCSAFTFamily; components; parameters::PCSAFTParams end
struct sPCSAFT <: sPCSAFTFamily; components; parameters::sPCSAFTParams end
struct ogSAFT <: ogSAFTFamily; components; parameters::ogSAFTParams end
