struct PCSAFTParams
    segment::Dict
    sigma::Dict
    epsilon::Dict
    epsilon_assoc::Dict
    bond_vol::Dict
    k::Dict
end

struct sPCSAFTParams
    segment::Dict
    sigma::Dict
    epsilon::Dict
    k::Dict
end

struct SAFTVRMieParams
    segment::Dict
    sigma::Dict
    epsilon::Dict
    lambdaA::Dict
    lambdaR::Dict
end

struct ogSAFTParams
    segment::Dict
    sigma::Dict
    epsilon::Dict
end
