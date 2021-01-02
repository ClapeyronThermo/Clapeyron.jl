abstract type Params end

struct PCSAFTParams <: Params
    segment::Dict #type: Dict{String,Float64}, can be changed to vector
    sigma::Dict #matrix of values
    epsilon::Dict #matrix of values
    epsilon_assoc::Dict #look on how to port this
    bond_vol::Dict #look on how to port this Dict{Set{Tuple{Set{String},String}},Float64}
    n_sites::Dict #dict of dicts
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

struct CKSAFTParams <: Params
    segment::Dict
    sigma::Dict
    epsilon::Dict
    c::Dict
    epsilon_assoc::Dict
    bond_vol::Dict
    n_sites::Dict
end

struct sCKSAFTParams <: Params
    segment::Dict
    sigma::Dict
    epsilon::Dict
    epsilon_assoc::Dict
    bond_vol::Dict
    n_sites::Dict
end

struct BACKSAFTParams <: Params
    segment::Dict
    sigma::Dict
    epsilon::Dict
    c::Dict
    alpha::Dict
end

struct SAFTVRSWParams <: Params
    segment::Dict
    sigma::Dict
    epsilon::Dict
    lambda::Dict
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

struct softSAFTParams <: Params
    segment::Dict #type: Dict{String,Float64}, can be changed to vector
    sigma::Dict #matrix of values
    epsilon::Dict #matrix of values
    epsilon_assoc::Dict #look on how to port this
    bond_vol::Dict #look on how to port this Dict{Set{Tuple{Set{String},String}},Float64}
    n_sites::Dict #dict of dicts
end

struct LJSAFTParams <: Params
    segment::Dict #type: Dict{String,Float64}, can be changed to vector
    b::Dict #matrix of values
    T::Dict #matrix of values
    epsilon_assoc::Dict #look on how to port this
    bond_vol::Dict #look on how to port this Dict{Set{Tuple{Set{String},String}},Float64}
    n_sites::Dict #dict of dicts
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
