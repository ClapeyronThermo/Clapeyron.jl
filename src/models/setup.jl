struct ModelOptions{P,S} 
    name::String
    parent::P
    siblings::S
    userlocations::Vector{String}
    mappings::Vector{Tuple{Union{String,Vector{String}},Symbol,Function}}
    has_components::Bool
    has_sites::Bool
    has_groups::Bool
    references::Vector{String}
    function ModelOptions(
            name::String;
            parent::P = nothing,
            siblings::S = Dict{String,ModelOptions}(),
            userlocations::Vector{String} = Vector{String}(),
            mappings::Vector{Tuple{Union{String,Vector{String}},Symbol,Function}} = Vector{Tuple{Union{String,Vector{String}},Symbol,Function}}(),
            has_components::Bool = true,
            has_sites::Bool = false,
            has_groups::Bool = false,
            references::Vector{String} = Vector{String}()
        ) where {P <: Union{ModelOptions,Nothing}, S <: Dict{String,ModelOptions}}
        return new{P,S}(
            name,
            parent,
            siblings,
            userlocations,
            mappings,
            has_components,
            has_sites,
            has_groups,
            references
        )
    end
end
