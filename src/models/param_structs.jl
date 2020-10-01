abstract type Params end

struct PCSAFTParams <: Params
    segment::Dict
    sigma::Dict
    epsilon::Dict
    epsilon_assoc::Dict
    bond_vol::Dict
    n_sites::Dict
end

struct sPCSAFTParams <: Params
    segment::Dict
    sigma::Dict
    epsilon::Dict
    epsilon_assoc::Dict
    bond_vol::Dict
    n_sites::Dict
    k::Dict
end

struct SAFTVRMieParams <: Params
    segment::Dict
    sigma::Dict
    epsilon::Dict
    lambdaA::Dict
    lambdaR::Dict
    epsilon_assoc::Dict
    bond_vol::Dict
    n_sites::Dict
end

struct ogSAFTParams <: Params
    segment::Dict
    sigma::Dict
    epsilon::Dict
    epsilon_assoc::Dict
    bond_vol::Dict
    n_sites::Dict
end
