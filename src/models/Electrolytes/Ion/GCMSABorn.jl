abstract type GCMSABornModel <: IonModel end

struct GCMSABornParam <: EoSParam
    shapefactor::SingleParam{Float64}
    segment::SingleParam{Float64}
    sigma::SingleParam{Float64}
    gc_sigma::SingleParam{Float64}
    sigma_born::SingleParam{Float64}
    gc_sigma_born::SingleParam{Float64}
    charge::SingleParam{Float64}
end

struct GCMSABorn{ϵ} <: GCMSABornModel
    components::Array{String,1}
    groups::GroupParam
    icomponents::UnitRange{Int}
    params::GCMSABornParam
    RSPmodel::ϵ
    references::Array{String,1}
end

export GCMSABorn

"""
    GCMSABorn(solvents::Array{String,1}, 
         ions::Array{String,1}; 
         RSPmodel=ConstW, 
         SAFTlocations=String[], 
         userlocations=String[], 
         verbose=false)

## Input parameters
- `sigma`: Single Parameter (`Float64`) - Hard-sphere diameter `[m]`
- `sigma_born`: Single Parameter (`Float64`) - Born Diameter `[m]`
- `charge`: Single Parameter (`Float64`) - Charge `[-]`

## Input models
- `RSPmodel`: Relative Static Permittivity Model

## Description
This function is used to create a group-contribution Mean Spherical Approximation-Born model used in SAFT-gamma E Mie
"""
function GCMSABorn(solvents,ions; RSPmodel=ConstW, userlocations=String[], verbose=false)
    groups = GroupParam(cat(solvents,ions,dims=1), ["SAFT/SAFTgammaMie/SAFTgammaMie_groups.csv"]; verbose=verbose)
    params = getparams(groups, ["SAFT/SAFTgammaMie/SAFTgammaMie_like.csv","SAFT/SAFTgammaMie/SAFTgammaMieE/","properties/molarmass_groups.csv"]; userlocations=userlocations,return_sites=false,ignore_missing_singleparams=["sigma_born","charge"], verbose=verbose)
    components = groups.components
    icomponents = 1:length(components)

    segment = params["vst"]
    shapefactor = params["S"]

    mix_segment!(groups,shapefactor.values,segment.values)

    sigma = params["sigma"]
    sigma.values .*= 1E-10
    gc_sigma = deepcopy(sigma)
    gc_sigma.values .^= 3
    gc_sigma.values .*= shapefactor.values .* segment.values
    gc_sigma.values .= cbrt.(gc_sigma.values)

    sigma_born = params["sigma_born"]
    sigma_born.values .*= 1E-10
    sigma_born.values[sigma_born.ismissingvalues] .= 1.07.*sigma.values[sigma_born.ismissingvalues]
    gc_sigma_born = deepcopy(sigma_born)
    gc_sigma_born.values .^= 3
    gc_sigma_born.values .*= shapefactor.values .* segment.values
    gc_sigma_born.values .= cbrt.(gc_sigma_born.values)

    charge = params["charge"]
    
    packagedparams = GCMSABornParam(shapefactor,segment,sigma,gc_sigma,sigma_born,gc_sigma_born,charge)

    references = String[]
    init_RSPmodel = RSPmodel(solvents,ions)


    model = GCMSABorn(components, groups, icomponents, packagedparams, init_RSPmodel,references)
    return model
end

function recombine_impl!(model::GCMSABornModel)
    groups = model.groups
    components = model.components

    sigma = deepcopy(model.params.sigma)
    segment = model.params.segment
    shapefactor = model.params.shapefactor

    sigma.values .^= 3
    sigma.values .*= shapefactor.values .* segment.values
    sigma.values .= cbrt.(sigma.values)

    model.params.gc_sigma.values .= sigma.values

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

function data(model::GCMSABornModel, V, T, z)
    ngroups = length(model.groups.flattenedgroups)
    v = model.groups.n_flattenedgroups
    Σz = sum(z)
    zg = zeros(eltype(sum(z)),ngroups)
    for k in 1:ngroups
        zg[k] = sum([z[i]*v[i][k] for i ∈ 1:length(model.groups.components)])
    end
    ng = sum(zg)
    return (zg, ng), dielectric_constant(model.RSPmodel, V, T, z)
end

function a_res(model::GCMSABornModel, V, T, z, _data=@f(data))
    return a_ion(model,V,T,z,_data)+a_born(model,V,T,z,_data)
end

function a_ion(model::GCMSABornModel, V, T, z, _data=@f(data))
    (zg, ∑zg), ϵ_r = _data
    σ = model.params.gc_sigma.values
    Z = model.params.charge.values
    igroups = 1:length(model.groups.flattenedgroups)
    iions = igroups[Z.!=0]
    if length(iions) == 0
        return zero(T+first(z))
    end
    ρg = N_A*sum(zg)/V
    Γ = @f(screening_length, ϵ_r, (zg, ∑zg))
    Δ = 1-π*ρg/6*sum(zg[i]*σ[i]^3 for i ∈ iions)/∑zg
    Ω = 1+π*ρg/(2*Δ)*sum(zg[i]*σ[i]^3/(1+Γ*σ[i]) for i ∈ iions)/∑zg
    Pn = ρg/Ω*sum(zg[i]*σ[i]*Z[i]/(1+Γ*σ[i]) for i ∈ iions)/∑zg

    U_GCMSA = -e_c^2*V/(4π*ϵ_0*ϵ_r)*(Γ*ρg*sum(zg[i]*Z[i]^2/(1+Γ*σ[i]) for i ∈ iions)/∑zg + π/(2Δ)*Ω*Pn^2)
    return (U_GCMSA+Γ^3*k_B*T*V/(3π))/(N_A*k_B*T*sum(z))
end

function screening_length(model::GCMSABornModel,V,T,z, ϵ_r=@f(data),zgdata = @f(data_msa))
    zg, ∑zg = zgdata
    ngroups = length(zg)

    σ = model.params.gc_sigma.values
    Z = model.params.charge.values
    igroups = 1:length(model.groups.flattenedgroups)
    iions = igroups[Z.!=0]

    #x = z ./ sum(z)
    ρg = N_A*∑zg/V
    Δ = 1-π*ρg/6*sum(zg[i]*σ[i]^3 for i ∈ iions)/∑zg

    Γold = (4π*e_c^2/(4π*ϵ_0*ϵ_r*k_B*T)*ρg*sum(zg[i]*Z[i]^2 for i ∈ iions)/∑zg)^(1/2)
    _0 = zero(Γold)
    Γnew = _0
    tol  = one(_0)
    iter = 1
    while tol>1e-12 && iter < 100
        Ω = 1+π*ρg/(2*Δ)*sum(zg[i]*σ[i]^3/(1+Γold*σ[i]) for i ∈ iions)/∑zg
        Pn = ρg/Ω*sum(zg[i]*σ[i]*Z[i]/(1+Γold*σ[i]) for i ∈ iions)/∑zg
        #Q = @. (Z-σ^2*Pn*(π/(2Δ)))./(1+Γold*σ)
        ∑Q2x = _0
        for i ∈ iions
            Qi = (Z[i]-σ[i]^2*Pn*(π/(2Δ)))/(1+Γold*σ[i])
            ∑Q2x += zg[i]*Qi^2
        end
        Γnew = sqrt(π*e_c^2*ρg/(4π*ϵ_0*ϵ_r*k_B*T)*∑Q2x/∑zg)
        tol = abs(1-Γnew/Γold)
        Γold = Γnew
        iter += 1
    end
    return Γnew
end

function a_born(model::GCMSABornModel, V, T, z,_data=@f(data)) 
    (zg, ∑zg), ϵ_r = _data  
    σ_born = model.params.sigma_born.values
    Z = model.params.charge.values
    igroups = 1:length(model.groups.flattenedgroups)
    iions = igroups[Z.!=0]
    if length(iions) == 0
        return zero(T+first(z))
    end
    

    return -e_c^2/(4π*ϵ_0*k_B*T*∑zg)*(1-1/ϵ_r)*sum(zg[i]*Z[i]^2/σ_born[i] for i ∈ iions)
end