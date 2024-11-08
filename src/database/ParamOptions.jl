const DEFAULT_N_SITES = Dict{String,String}()

const IGNORE_HEADERS = ["dipprnumber", "smiles", "cas"]

"""
    ParamOptions(;kwargs...)
Struct containing all the options related to parameter parsing:
* `userlocations = String[]`: List of used-defined locations to search.
* `group_userlocations = String[]`: List of used-defined group locations to search.
* `asymmetricparams::Vector{String} = String[]`: List of pair parameters that follow that `param[i,j] ≠ param[j,i]`. if not set on asymmetric pairs, the asymmetric values will be overwritten!
* `ignore_headers::Vector{String} = ["dipprnumber", "smiles"]`: List of ignored headers.
* `ignore_missing_singleparams::Vector{String} = String[]`: List of parameters where checking for missing parameter values (in `SingleParam`) or the diagonal (on `PairParam`) are ignored.
* `verbose::Bool = false`: If `true`, show all operations done by `getparams` displayed in the terminal. this includes the warnings emmited by `CSV.jl`
* `species_columnreference::String ="species"`: column name to check for components. in pair and association params, it will check for `#species#1` and `#species#2`, where `#species#` is the value of this option.
* `site_columnreference::String ="site"`: column name to check for sites in association params, it will check for `#site#1` and `#site#2`, where `#site#` is the value of this option.
* `normalisecomponents::Bool = true`: If `true`, performs normalization of strings, on the CSV and input components
* `n_sites_columns::Dict{String,String} = Dict( "e" => "n_e","e1" => "n_e1","e2" => "n_e2","H" => "n_H")`: dictionary to look number of sites. the number of sites is stored as columns in a single parameter csv file. for example, the number of sites of name `e` will be looked on the column `n_e`
* `return_sites::Bool = true`: If set to false, association params will be ignored and sites will not be created, even if they exist in the list of locations.
* `component_delimiter::String = "~|~"`: When there are multiple component names to match, seperate them by this delimiter.
"""
Base.@kwdef struct ParamOptions{𝕌,𝔾}
    userlocations::𝕌 = String[]
    group_userlocations::𝔾 = String[]
    asymmetricparams::Vector{String}= String[]
    ignore_missing_singleparams::Vector{String} = String[]
    ignore_headers::Vector{String} = IGNORE_HEADERS
    verbose::Bool = false
    species_columnreference::String ="species"
    source_columnreference::String = "source"
    site_columnreference::String = "site"
    normalisecomponents::Bool = true
    n_sites_columns::Dict{String,String} = DEFAULT_N_SITES
    return_sites::Bool = true
    component_delimiter::String = "~|~"
end

const DefaultOptions = ParamOptions()
const DefaultGroupOptions = ParamOptions(ignore_missing_singleparams = ["intragroups"])

export ParamOptions
