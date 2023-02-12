abstract type GCMSAModel <: IonModel end

struct GCMSAParam <: EoSParam
    shapefactor::SingleParam{Float64}
    segment::SingleParam{Float64}
    sigma::SingleParam{Float64}
    gc_sigma::SingleParam{Float64}
    charge::SingleParam{Float64}
end

struct GCMSA{ϵ} <: GCMSAModel
    components::Array{String,1}
    groups::GroupParam
    icomponents::UnitRange{Int}
    isolvents::UnitRange{Int}
    iions::UnitRange{Int}
    params::GCMSAParam
    RSPmodel::ϵ
    absolutetolerance::Float64
    references::Array{String,1}
end

export GCMSA
function GCMSA(solvents,salts,ions; RSPmodel=ConstW, SAFTlocations=String[], userlocations=String[], ideal_userlocations=String[], verbose=false)
    groups = GroupParam(cat(solvents,ions,dims=1), ["SAFT/SAFTgammaMie/SAFTgammaMie_groups.csv"]; verbose=verbose)
    params = getparams(groups, ["SAFT/SAFTgammaMie/SAFTgammaMie_like.csv","SAFT/SAFTgammaMie/SAFTgammaMieE/","properties/molarmass_groups.csv"]; userlocations=userlocations,return_sites=false,ignore_missing_singleparams=["sigma_born","charge"], verbose=verbose)
    components = groups.components

    segment = params["vst"]
    shapefactor = params["S"]

    mix_segment!(groups,shapefactor.values,segment.values)

    sigma = params["sigma"]
    sigma.values .*= 1E-10
    gc_sigma = deepcopy(sigma)
    gc_sigma.values .^= 3
    gc_sigma.values .*= shapefactor.values .* segment.values
    gc_sigma.values .= cbrt.(gc_sigma.values)

    charge = params["charge"]

    components = groups.components
    icomponents = 1:length(components)
    isolvents = 1:length(solvents)
    iions = (length(solvents)+1):length(components)
    
    packagedparams = GCMSAParam(shapefactor,segment,sigma,gc_sigma,charge)

    references = String[]
    if RSPmodel !== nothing
        init_RSPmodel = RSPmodel(solvents,salts)
    else
        init_RSPmodel = nothing
    end

    model = GCMSA(components, groups, icomponents, isolvents, iions, packagedparams, init_RSPmodel, 1e-12,references)
    return model
end

function recombine_impl!(model::GCMSAModel)
    groups = model.groups
    components = model.components

    sigma = deepcopy(model.params.sigma)
    segment = model.params.segment
    shapefactor = model.params.shapefactor

    sigma.values .^= 3
    sigma.values .*= shapefactor.values .* segment.values
    sigma.values .= cbrt.(sigma.values)

    model.params.gc_sigma.values .= sigma.values
    return model
end

function data(model::GCMSAModel, V, T, z)
    return @f(data_msa),@f(data_rsp)
end

function data_msa(model::GCMSAModel, V, T, z)
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

function data_rsp(model::GCMSAModel, V, T, z)
    return dielectric_constant(model.RSPmodel, V, T, z)
end

function a_res(model::GCMSAModel, V, T, z, _data=@f(data))
    (zg, ∑zg), ϵ_r = _data
    ngroups = length(zg)
    if ngroups == 0
        return zero(V+T+first(z))
    end
    σ = model.params.gc_sigma.values
    Z = model.params.charge.values

    ρg = N_A*sum(zg)/V
    Γ = @f(screening_length, _data)
    Δ = 1-π*ρg/6*sum(zg[i]*σ[i]^3 for i ∈ 1:ngroups)/∑zg
    Ω = 1+π*ρg/(2*Δ)*sum(zg[i]*σ[i]^3/(1+Γ*σ[i]) for i ∈ 1:ngroups)/∑zg
    Pn = ρg/Ω*sum(zg[i]*σ[i]*Z[i]/(1+Γ*σ[i]) for i ∈ 1:ngroups)/∑zg

    U_GCMSA = -e_c^2*V/(4π*ϵ_0*ϵ_r)*(Γ*ρg*sum(zg[i]*Z[i]^2/(1+Γ*σ[i]) for i ∈ 1:ngroups)/∑zg + π/(2Δ)*Ω*Pn^2)
    return (U_GCMSA+Γ^3*k_B*T*V/(3π))/(N_A*k_B*T*sum(z))
end

function screening_length(model::GCMSAModel,V,T,z, _data=@f(data))
    (zg, ∑zg), ϵ_r = _data
    ngroups = length(zg)

    σ = model.params.gc_sigma.values
    Z = model.params.charge.values
    #x = z ./ sum(z)
    ρg = N_A*∑zg/V
    Δ = 1-π*ρg/6*sum(zg[i]*σ[i]^3 for i ∈ 1:ngroups)/∑zg

    Γold = (4π*e_c^2/(4π*ϵ_0*ϵ_r*k_B*T)*ρg*sum(zg[i]*Z[i]^2 for i ∈ 1:ngroups)/∑zg)^(1/2)
    _0 = zero(Γold)
    Γnew = _0
    tol  = one(_0)
    iter = 1
    while tol>1e-8 && iter < 100
        Ω = 1+π*ρg/(2*Δ)*sum(zg[i]*σ[i]^3/(1+Γold*σ[i]) for i ∈ 1:ngroups)/∑zg
        Pn = ρg/Ω*sum(zg[i]*σ[i]*Z[i]/(1+Γold*σ[i]) for i ∈ 1:ngroups)/∑zg
        #Q = @. (Z-σ^2*Pn*(π/(2Δ)))./(1+Γold*σ)
        ∑Q2x = _0
        for i ∈ 1:ngroups
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