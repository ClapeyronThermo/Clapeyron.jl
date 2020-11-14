function system(components::Array{String,1}, method::String; kwargs...)
    # possible kwargs... are filepaths for
    # customdatabase_like, customdatabase_unlike, customdatabase_assoc
    set_components = [Set([i]) for i in components]
    if method == "PCSAFT"
        raw_params = retrieveparams(components, method; kwargs...)
        params = create_PCSAFTParams(raw_params)
        sites = extractsites(params.n_sites)
        model = PCSAFT(set_components, sites, params)
    elseif method == "SAFTVRMie"
        raw_params = retrieveparams(components, method; kwargs...)
        model = SAFTVRMie(set_components, create_SAFTVRMieParams(raw_params))
    elseif method == "SAFTVRQMie"
        raw_params = retrieveparams(components, method; kwargs...)
        model = SAFTVRQMie(set_components, create_SAFTVRQMieParams(raw_params))
    elseif method == "ogSAFT"
        raw_params = retrieveparams(components, method; kwargs...)
        params = create_ogSAFTParams(raw_params)
        model = ogSAFT(set_components, params)
    elseif method == "sPCSAFT"
        raw_params = retrieveparams(components, method; kwargs...)
        params = create_sPCSAFTParams(raw_params)
        model = sPCSAFT(set_components, params)
    elseif method == "vdW"
        raw_params = retrieveparams(components, method; kwargs...)
        params = create_vdWParams(raw_params)
        model = vdW(set_components, params)
    elseif method == "RK"
        raw_params = retrieveparams(components, method; kwargs...)
        params = create_RKParams(raw_params)
        model = RK(set_components, params)
    elseif method == "SRK"
        raw_params = retrieveparams(components, method; kwargs...)
        params = create_SRKParams(raw_params)
        model = SRK(set_components, params)
    elseif method == "PR"
        raw_params = retrieveparams(components, method; kwargs...)
        params = create_PRParams(raw_params)
        model = PR(set_components, params)
    elseif method == "CPA"
        raw_params = retrieveparams(components, method; kwargs...)
        params = create_CPAParams(raw_params)
        model = CPA(set_components, params)
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
        params = create_SAFTgammaMie(raw_params)
        model = SAFTgammaMie(components, groups, group_multiplicities, params)
    else
        error("Method definition incorrect.")
    end
    return model
end
