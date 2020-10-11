abstract type EoS end
abstract type SAFT <: EoS end
abstract type Ideal end

abstract type PCSAFTFamily <: SAFT end
abstract type sPCSAFTFamily <: SAFT end
abstract type ogSAFTFamily <: SAFT end
abstract type SAFTVRMieFamily <: SAFT end

struct SAFTVRMie <: SAFTVRMieFamily; components; params::SAFTVRMieParams end
struct PCSAFT <: PCSAFTFamily; components; params::PCSAFTParams end
struct sPCSAFT <: sPCSAFTFamily; components; params::sPCSAFTParams end
struct ogSAFT <: ogSAFTFamily; components; params::ogSAFTParams end

struct Monomer <: Ideal; components; params end
struct Reid <: Ideal; components; params end
