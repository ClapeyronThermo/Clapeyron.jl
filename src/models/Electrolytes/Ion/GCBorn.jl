abstract type GCBornModel <: EoSModel end

struct GCBornParam <: EoSParam
    shapefactor::SingleParam{Float64}
    segment::SingleParam{Float64}
    sigma::SingleParam{Float64}
    sigma_born::SingleParam{Float64}
    gc_sigma_born::SingleParam{Float64}
    charge::SingleParam{Float64}
end
struct GCBorn{ϵ} <: GCBornModel
    components::Array{String,1}
    groups::GroupParam
    isolvents::UnitRange{Int}
    iions::UnitRange{Int}
    params::GCBornParam
    RSPmodel::ϵ
    references::Array{String,1}
end

export GCBorn
function GCBorn(solvents,salts,ions; RSPmodel=ConstW, SAFTlocations=String[], userlocations=String[], ideal_userlocations=String[], verbose=false)
    groups = GroupParam(cat(solvents,ions,dims=1), ["SAFT/SAFTgammaMie/SAFTgammaMie_groups.csv"]; verbose=verbose)
    params = getparams(groups, ["SAFT/SAFTgammaMie/SAFTgammaMie_like.csv","SAFT/SAFTgammaMie/SAFTgammaMieE/","properties/molarmass_groups.csv"]; userlocations=userlocations,return_sites=false,ignore_missing_singleparams=["sigma_born","charge"], verbose=verbose)
    components = groups.components
    isolvents = 1:length(solvents)
    iions = (length(solvents)+1):length(components)

    segment = params["vst"]
    shapefactor = params["S"]
    sigma = params["sigma"]
    sigma.values .*= 1e-10

    mix_segment!(groups,shapefactor.values,segment.values)

    sigma_born = params["sigma_born"]
    sigma_born.values .*= 1E-10
    sigma_born.values[sigma_born.ismissingvalues] .= 1.07.*sigma.values[sigma_born.ismissingvalues]
    gc_sigma_born = deepcopy(sigma_born)
    gc_sigma_born.values .^= 3
    gc_sigma_born.values .*= shapefactor.values .* segment.values
    gc_sigma_born.values .= cbrt.(gc_sigma_born.values)

    charge = params["charge"]
    
    packagedparams = GCBornParam(shapefactor,segment,sigma,sigma_born,gc_sigma_born,charge)

    references = String[]
    if RSPmodel !== nothing
        init_RSPmodel = RSPmodel(solvents,salts)
    else
        init_RSPmodel = nothing
    end

    model = GCBorn(components, groups, isolvents, iions, packagedparams, init_RSPmodel, references)
    return model
end

function recombine_impl!(model::GCBornModel)
    groups = model.groups
    components = model.components

    model.params.sigma_born.values[model.params.sigma_born.ismissingvalues] .= 1.07.*model.params.sigma.values[model.params.sigma_born.ismissingvalues]

    sigma_born = deepcopy(model.params.sigma_born)
    segment = model.params.segment
    shapefactor = model.params.shapefactor

    sigma_born.values .^= 3
    sigma_born.values .*= shapefactor.values .* segment.values
    sigma_born.values .= cbrt.(sigma_born.values)

    model.params.gc_sigma_born.values .= sigma_born.values
    return model
end

function data(model::GCBornModel, V, T, z)
    return @f(data_rsp)
end

function data_born(model::GCBornModel, V, T, z)
    groups = model.groups
    ngroups = length(groups.flattenedgroups)
    v = groups.n_flattenedgroups
    zg = zeros(eltype(sum(z)),ngroups)
    for k in 1:ngroups
        zg[k] = sum([z[i]*v[i][k] for i ∈ 1:length(groups.components)])
    end
    ng = sum(zg)
    return  zg, ng
end

function data_rsp(model::GCBornModel, V, T, z)
    return dielectric_constant(model.RSPmodel, V, T, z)
end


function a_res(model::GCBornModel, V, T, z,_data=@f(data))   
    ngroups = length(model.groups.flattenedgroups)
    if ngroups == 0
        return zero(T+first(z))
    end
    σ_born = model.params.sigma_born.values
    Z = model.params.charge.values
    ϵ_r = _data
    zg, ng = data_born(model,V,T,z)

    return -e_c^2/(4π*ϵ_0*k_B*T*sum(z))*(1-1/ϵ_r)*sum(zg[i]*Z[i]^2/σ_born[i] for i ∈ 1:ngroups)
end