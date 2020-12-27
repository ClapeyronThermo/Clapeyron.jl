function system(components::Array{String,1}, method::String, ideal="Basic"; kwargs...)
    # possible kwargs... are filepaths for
    # customdatabase_like, customdatabase_unlike, customdatabase_assoc
    set_components = [Set([i]) for i in components]
    if method == "PCSAFT"
        raw_params = retrieveparams(components, method, ideal; kwargs...)
        params = create_PCSAFTParams(raw_params)
        sites = extractsites(params.n_sites)
        ideal_model = create_IdealParams(set_components, raw_params,ideal)
        model = PCSAFT(set_components, sites, params, ideal_model)
    elseif method == "SAFTVRMie"
        raw_params = retrieveparams(components, method, ideal; kwargs...)
        ideal_model = create_IdealParams(set_components, raw_params,ideal)
        model = SAFTVRMie(set_components, create_SAFTVRMieParams(raw_params),ideal_model)
    elseif method == "SAFTVRQMie"
        raw_params = retrieveparams(components, method, ideal; kwargs...)
        ideal_model = create_IdealParams(set_components, raw_params,ideal)
        model = SAFTVRQMie(set_components, create_SAFTVRQMieParams(raw_params),ideal_model)
    elseif method == "SAFTVRSW"
        raw_params = retrieveparams(components, method, ideal; kwargs...)
        ideal_model = create_IdealParams(set_components, raw_params,ideal)
        model = SAFTVRSW(set_components, create_SAFTVRSWParams(raw_params),ideal_model)
    elseif method == "CKSAFT"
        raw_params = retrieveparams(components, method, ideal; kwargs...)
        ideal_model = create_IdealParams(set_components, raw_params,ideal)
        model = CKSAFT(set_components, create_CKSAFTParams(raw_params),ideal_model)
    elseif method == "ogSAFT"
        raw_params = retrieveparams(components, method, ideal; kwargs...)
        params = create_ogSAFTParams(raw_params)
        ideal_model = create_IdealParams(set_components, raw_params,ideal)
        model = ogSAFT(set_components, params,ideal_model)
    elseif method == "softSAFT"
        raw_params = retrieveparams(components, method, ideal; kwargs...)
        params = create_softSAFTParams(raw_params)
        ideal_model = create_IdealParams(set_components, raw_params,ideal)
        model = softSAFT(set_components, params,ideal_model)
    elseif method == "sPCSAFT"
        raw_params = retrieveparams(components, method, ideal; kwargs...)
        ideal_model = create_IdealParams(set_components, raw_params,ideal)
        params = create_sPCSAFTParams(raw_params)
        model = sPCSAFT(set_components, params,ideal_model)
    elseif method == "vdW"
        raw_params = retrieveparams(components, method, ideal; kwargs...)
        params = create_vdWParams(raw_params)
        ideal_model = create_IdealParams(set_components, raw_params,ideal)
        model = vdW(set_components, params, ideal_model)
    elseif method == "RK"
        raw_params = retrieveparams(components, method, ideal; kwargs...)
        ideal_model = create_IdealParams(set_components, raw_params,ideal)
        params = create_RKParams(raw_params)
        model = RK(set_components, params, ideal_model)
    elseif method == "SRK"
        raw_params = retrieveparams(components, method, ideal; kwargs...)
        ideal_model = create_IdealParams(set_components, raw_params,ideal)
        params = create_SRKParams(raw_params)
        model = SRK(set_components, params, ideal_model)
    elseif method == "PR"
        raw_params = retrieveparams(components, method, ideal; kwargs...)
        ideal_model = create_IdealParams(set_components, raw_params,ideal)
        params = create_PRParams(raw_params)
        model = PR(set_components, params, ideal_model)
    elseif method == "CPA"
        raw_params = retrieveparams(components, method, ideal; kwargs...)
        ideal_model = create_IdealParams(set_components, raw_params,ideal)
        params = create_CPAParams(raw_params)
        model = CPA(set_components, params, ideal_model)
    else
        error("Method definition incorrect.")
    end
    return model
end

function system(component::String, method::String,ideal="Basic"; kwargs...)
    return system([component], method; kwargs...)
end

function system(group_multiplicities::Dict, method::String, ideal="Basic"; kwargs...)
    # possible kwargs... are filepaths for
    # customdatabase_like, customdatabase_unlike, customdatabase_assoc
    components = collect(keys(group_multiplicities))
    set_components = [i for i in components]
    groups = union(vcat([collect(keys(i)) for i in values(group_multiplicities)]...))
    string_groups = [collect(j)[1] for j in groups]
    if method == "SAFTgammaMie"
        raw_params = retrieveparams(string_groups, method, ideal; kwargs...)
        ideal_model = create_IdealParams(set_components, raw_params,ideal)
        params = create_SAFTgammaMie(raw_params)
        model = SAFTgammaMie(components, groups, group_multiplicities, params, ideal_model)
    else
        error("Method definition incorrect.")
    end
    return model
end
