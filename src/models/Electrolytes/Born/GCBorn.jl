abstract type GCBornModel <: EoSModel end

struct GCBorn{ϵ} <: GCBornModel
    components::Array{String,1}
    groups::GroupParam
    icomponents::UnitRange{Int}
    isolvents::UnitRange{Int}
    iions::UnitRange{Int}
    params::BornParam
    RSPmodel::ϵ
    absolutetolerance::Float64
    references::Array{String,1}
end

@registermodel GCBorn
export GCBorn
function GCBorn(solvents,salts,ions; RSPmodel=ConstW, SAFTlocations=String[], userlocations=String[], ideal_userlocations=String[], verbose=false)
    groups = GroupParam(cat(solvents,ions,dims=1), ["SAFT/SAFTgammaMie/SAFTgammaMie_groups.csv"]; verbose=verbose)
    params = getparams(groups, ["SAFT/SAFTgammaMie/SAFTgammaMie_like.csv","SAFT/SAFTgammaMie/SAFTgammaMieE/","properties/molarmass_groups.csv"]; userlocations=userlocations, verbose=verbose)
    components = groups.components

    gc_segment = params["vst"]
    shapefactor = params["S"]

    mix_segment!(groups,shapefactor.values,gc_segment.values)

    gc_sigma_born = params["sigma_born"]
    gc_sigma_born.values .*= 1E-10
    gc_sigma_born.values .^= 3
    gc_sigma_born.values .*= shapefactor.values .* gc_segment.values
    gc_sigma_born.values .= cbrt.(gc_sigma_born.values)

    charge = params["charge"]

    components = groups.components
    icomponents = 1:length(components)
    isolvents = 1:length(solvents)
    iions = (length(solvents)+1):length(components)
    
    packagedparams = BornParam(gc_sigma_born,charge)

    references = String[]
    if RSPmodel !== nothing
        init_RSPmodel = RSPmodel(solvents,salts)
    else
        init_RSPmodel = nothing
    end

    model = GCBorn(components, groups, icomponents, isolvents, iions, packagedparams, init_RSPmodel, 1e-12,references)
    return model
end

function data(model::GCBornModel, V, T, z)
    return @f(data_born),@f(data_rsp)
end

function data_born(model::GCBornModel, V, T, z)
    ngroups = length(model.groups.flattenedgroups)
    v = model.groups.n_flattenedgroups
    Σz = sum(z)
    zg = zeros(eltype(sum(z)),ngroups)
    for k in 1:ngroups
        zg[k] = sum([z[i]*v[i][k] for i ∈ model.icomponents])
    end
    ng = sum(zg)
    return  zg, ng
end

function data_rsp(model::GCBornModel, V, T, z)
    return dielectric_constant(model.rspmodel, V, T, z)
end


function a_res(model::GCBornModel, V, T, z,_data=@f(data))
    σ_born = model.params.sigma_born.values
    Z = model.params.charge.values
    (zg,ng),ϵ_r = _data
    ngroups = length(zg)

    return -e_c^2/(4π*ϵ_0*k_B*T*sum(z))*(1-1/ϵ_r)*sum(zg[i]*Z[i]^2/σ_born[i] for i ∈ 1:ngroups)
end