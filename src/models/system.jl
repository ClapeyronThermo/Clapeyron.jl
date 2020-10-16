function system(components::Array{String,1}, method::String; kwargs...)
    # possible kwargs... are filepaths for
    # customdatabase_like, customdatabase_unlike, customdatabase_assoc
    set_components = [Set([i]) for i in components]
    if method == "PCSAFT"
        raw_params = retrieveparams(components, method; kwargs...)
        model = PCSAFT(set_components, create_PCSAFTParams(raw_params))
    elseif method == "SAFTVRMie"
        raw_params = retrieveparams(components, method; kwargs...)
        model = SAFTVRMie(set_components, create_SAFTVRMieParams(raw_params))
    elseif method == "ogSAFT"
        raw_params = retrieveparams(components, method; kwargs...)
        model = ogSAFT(set_components, create_ogSAFTParams(raw_params))
    elseif method == "sPCSAFT"
        raw_params = retrieveparams(components, method; kwargs...)
        model = sPCSAFT(set_components, create_sPCSAFTParams(raw_params))
    else
        error("Method definition incorrect.")
    end
    return model
end

function system(component::String, method::String; kwargs...)
    return system([component], method; kwargs...)
end

function system(group_multiplicities::Dict, method::String; kwargs...)
    # possible kwargs... are filepaths for
    # customdatabase_like, customdatabase_unlike, customdatabase_assoc
    components = collect(keys(group_multiplicities))
    groups = union(vcat([collect(keys(i)) for i in values(group_multiplicities)]...))
    string_groups = [collect(j)[1] for j in groups]
    if method == "SAFTgammaMie"
        raw_params = retrieveparams(string_groups, method; kwargs...)
        model = SAFTgammaMie(components, groups, group_multiplicities, create_SAFTgammaMie(raw_params))
    else
        error("Method definition incorrect.")
    end
    return model
end
