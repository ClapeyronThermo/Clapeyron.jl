const DEFAULT_N_SITES = Dict{String,String}(
    "e" => "n_e",
    "e1" => "n_e1",
    "e2" => "n_e2",
    "H" => "n_H"
)

Base.@kwdef struct ParamOptions
    userlocations::Vector{String} =String[]
    usergrouplocations::Vector{String} = String[]
    asymmetricparams::Vector{String}=String[]
    ignore_missing_singleparams::Bool=false 
    verbose::Bool=false
    species_columnreference::String="species"
    source_columnreference::String="source"
    site_columnreference::String="site"
    group_columnreference::String="groups" 
    normalisecomponents::Bool=true
    n_sites_columns::Dict{String,String} = DEFAULT_N_SITES
    return_sites::Bool = true
end