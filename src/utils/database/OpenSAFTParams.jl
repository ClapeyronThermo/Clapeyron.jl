abstract type OpenSAFTParams end

struct SingleParams{T} <: OpenSAFTParams
    name::String
    values::Array{T,1}
    components::Array{String,1}
    model::String
    sources::Array{String,1}
end

struct PairParams{T} <: OpenSAFTParams
    name::String
    values::Array{T,2}
    components::Array{String,1}
    model::String
    sources::Array{String,1}
end

struct AssocParams{T} <: OpenSAFTParams
    name::String
    values::Array{Array{T,2},2}
    components::Array{String,1}
    sites::Array{Array{String,1},1}
    model::String
    sources::Array{String,1}
end

struct GroupParams <: OpenSAFTParams
    components::Array{String,1}
    groups::Array{Array{String,1},1}
    groupmultiplicities::Array{Array{Int64,1},1}
    model::String
end
