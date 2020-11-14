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

struct SAFTVRQMieParams <: Params
    MolarMass::Dict
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

struct SAFTgammaMieParams <: Params
    segment::Dict
    shapefactor::Dict
    lambda_a::Dict
    lambda_r::Dict
    sigma::Dict
    epsilon::Dict
    epsilon_assoc::DefaultDict
    bond_vol::DefaultDict
    n_sites::DefaultDict
end

struct vdWParams <: Params
    a::Dict
    b::Dict
    Tc::Dict
end

struct RKParams <: Params
    a::Dict
    b::Dict
    Tc::Float64
end

struct SRKParams <: Params
    a::Dict
    b::Dict
    Tc::Dict
    acentric_fac::Dict
end

struct PRParams <: Params
    a::Dict
    b::Dict
    Tc::Dict
    acentric_fac::Dict
end

struct CPAParams <: Params
    a::Dict
    b::Dict
    c1::Dict
    Tc::Dict
    epsilon_assoc::Dict
    bond_vol::Dict
    n_sites::Dict
end
