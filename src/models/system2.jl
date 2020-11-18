#=
the idea is the following:
system is a thin wrapper to CreateModel, where CreateModel is dispatched to each EoS
If a string is provided,conversion can be done to a proper type.
#this conversion is done by the function eos_from_text, see in this file how to specify new eos names
=#

#you can pass a collection of strings here, no need to specify.
function System(components, method::Type{T}; kwargs...) where T <: EoS
    return CreateModel(method,components,kwargs...)
end

function System(components, method::Union{String,Symbol}; kwargs...)
    return CreateModel(eos_from_text(method),components,kwargs...)
end

#if only a string is passed,convert to a collection.
function System(components::String, method; kwargs...)
    return System([components], method; kwargs...)
end

eos_from_text(str::String) = eos_from_text(Symbol(str))
eos_from_text(sym::Symbol) = eos_from_text(Val{sym})

#catch all
function eos_from_text(val::Type{Val{T}}) where T 
    return error("Method definition incorrect.")
end

#actual dispaches, adding a new synonym to a equation of state is as symple as adding a new 
#dispatch here, without need to modify the system function

#this is a highly dynamic dispatch, it is slow in the context of high performance calculations (4 μs, does not inline)
#but, on the other part, allows crazy things like this
#good thing you only create a few models every session
eos_from_text(val::Type{Val{:PCSAFT}}) = PCSAFT
eos_from_text(val::Type{Val{:SAFTVRMie}}) = SAFTVRMie
eos_from_text(val::Type{Val{:SAFTVRQMie}}) = SAFTVRQMie
eos_from_text(val::Type{Val{:ogSAFT}}) = ogSAFT
eos_from_text(val::Type{Val{:sPCSAFT}}) = sPCSAFT
eos_from_text(val::Type{Val{:vdW}}) = vdW
eos_from_text(val::Type{Val{:VdW}}) = vdW
eos_from_text(val::Type{Val{:RK}}) = RK
eos_from_text(val::Type{Val{:SRK}}) = SRK
eos_from_text(val::Type{Val{:PR}}) = PR
eos_from_text(val::Type{Val{:CPA}}) = CPA
eos_from_text(val::Type{Val{:SAFTgammaMie}}) = SAFTgammaMie
eos_from_text(val::Type{Val{:SAFTγMie}}) = SAFTgammaMie


function CreateModel(method::Type{PCSAFT},components;kwargs...)
    set_components = [Set([i]) for i in components]
    raw_params = retrieveparams(components, "PCSAFT"; kwargs...)
    params = create_PCSAFTParams(raw_params)
    sites = extractsites(params.n_sites)
    model = PCSAFT(set_components, sites, params)
end

function CreateModel(method::Type{SAFTVRMie},components;kwargs...)
    set_components = [Set([i]) for i in components]
    raw_params = retrieveparams(components,"SAFTVRMie"; kwargs...)
    model = SAFTVRMie(set_components, create_SAFTVRMieParams(raw_params))
end

function CreateModel(method::Type{SAFTVRQMie},components;kwargs...)
    set_components = [Set([i]) for i in components]
    raw_params = retrieveparams(components, "SAFTVRQMie"; kwargs...)
    model = SAFTVRQMie(set_components, create_SAFTVRQMieParams(raw_params))
end

function CreateModel(method::Type{ogSAFT},components;kwargs...)
    set_components = [Set([i]) for i in components]
    raw_params = retrieveparams(components, "ogSAFT"; kwargs...)
    params = create_ogSAFTParams(raw_params)
    model = ogSAFT(set_components, params)
end

function CreateModel(method::Type{sPCSAFT},components;kwargs...)
    set_components = [Set([i]) for i in components]
    raw_params = retrieveparams(components, "sPCSAFT"; kwargs...)
    params = create_sPCSAFTParams(raw_params)
    model = sPCSAFT(set_components, params)
end

function CreateModel(method::Type{vdW},components;kwargs...)
    set_components = [Set([i]) for i in components]
    raw_params = retrieveparams(components, "vdW"; kwargs...)
    params = create_vdWParams(raw_params)
    model = vdW(set_components, params)
end

function CreateModel(method::Type{RK},components;kwargs...)
    set_components = [Set([i]) for i in components]
    raw_params = retrieveparams(components, "RK"; kwargs...)
    params = create_RKParams(raw_params)
    model = RK(set_components, params)
end

function CreateModel(method::Type{SRK},components;kwargs...)
    set_components = [Set([i]) for i in components]
    raw_params = retrieveparams(components,"SRK"; kwargs...)
    params = create_SRKParams(raw_params)
    model = SRK(set_components, params)
end

function CreateModel(method::Type{PR},components;kwargs...)
    set_components = [Set([i]) for i in components]
    raw_params = retrieveparams(components, "PR"; kwargs...)
    params = create_PRParams(raw_params)
    model = PR(set_components, params)
end

function CreateModel(method::Type{CPA},components;kwargs...)
    set_components = [Set([i]) for i in components]
    raw_params = retrieveparams(components, "CPA"; kwargs...)
    params = create_CPAParams(raw_params)
    model = CPA(set_components, params)
end

#here, instead of a collection of strings, a key=> value collection is passed (Dict,NamedTuple).
function CreateModel(method::Type{SAFTgammaMie},group_multiplicities;kwargs...)
    # possible kwargs... are filepaths for
    # customdatabase_like, customdatabase_unlike, customdatabase_assoc
    components = collect(keys(group_multiplicities))
    groups = union(vcat([collect(keys(i)) for i in values(group_multiplicities)]...))
    string_groups = [collect(j)[1] for j in groups]
    raw_params = retrieveparams(string_groups, "SAFTgammaMie"; kwargs...)
    params = create_SAFTgammaMie(raw_params)
    model = SAFTgammaMie(components, groups, group_multiplicities, params)
    return model
end

