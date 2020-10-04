function system(components::Array{String,1}, method::String; kwargs...)
    # possible kwargs... are filepaths for
    # customdatabase_like, customdatabase_unlike, customdatabase_assoc
    set_components = [Set([components[i]]) for i in 1:length(components)]
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
    end
    return model
end

function system(component::String, method::String; kwargs...)
    return system([component], method; kwargs...)
end
