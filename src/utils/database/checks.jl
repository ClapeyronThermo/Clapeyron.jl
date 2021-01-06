function checkfor_clashingheaders(filepaths::Array{String,1})
    # Raises an error if the header of any assoc parameter clashes with a non-assoc parameter
    headerparams = []
    headerparams_assoc = []
    for filepath in filepaths
        type = readtype(filepath)
        if type == singledata || type == pairdata
            append!(headerparams, readheader(filepath))
        elseif type == assocdata
            append!(headerparams_assoc, readheader(filepath))
        end
    end
    clashingheaders = intersect(headerparams, headerparams_assoc)
    !isempty(clashingheaders) && error("Headers ", clashingheaders, " appear in both loaded assoc and non-assoc files.")
end

function checkfor_singlecompleteness(filepaths::Array{String,1}, components::Array{String,1})
    # Raises an error if any component is not present in any single databases.
end
