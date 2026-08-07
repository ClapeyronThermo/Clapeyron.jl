abstract type MMModel <: RSPModel end

struct MM1Param <: EoSParam
    gamma::PairParam{Float64}
    theta::PairParam{Float64}
    coordz::PairParam{Float64}
    mu::SingleParam{Float64}
    polarizability::SingleParam{Float64}
end

struct MM1 <: MMModel
    components::Array{String,1}
    params::MM1Param
    references::Array{String,1}
end

IonDependency(model::MM1) = DependentIonModel(model)

function MM1(solvents,ions; userlocations = String[], verbose::Bool=false)
    components = vcat(solvents,ions)
    params = getparams(components, ["Electrolytes/RSP/MM1_like.csv","Electrolytes/RSP/MM1_unlike.csv"]; userlocations=userlocations, verbose=verbose)
    params["gamma"].values .*= pi/180
    gamma = params["gamma"]
    params["theta"].values .*= pi/180
    theta = params["theta"]
    coordz = params["coordz"]
    params["mu"].values .*= 1. /(299792458)*1e-21
    mu = params["mu"]
    params["polarizability"].values .*= 1e-40
    polarizability = params["polarizability"]
    packagedparams = MM1Param(PairParam(gamma),PairParam(theta),PairParam(coordz),mu,polarizability)
    references = String[]
    model = MM1(components, packagedparams, references)
    return model
end

export MM1

#implementation of a_res_minus_assoc for most SAFT/CPA models
function a_res_minus_assoc(model::PCSAFTModel,V,T,z,_data)
    @f(a_hc,_data) + @f(a_disp,_data)
end

function a_res_minus_assoc(model ::PCPSAFTModel, V, T, z,_data)
    return @f(a_hc,_data) + @f(a_disp,_data) + @f(a_polar,_data)
end

function a_res_minus_assoc(model ::QPCPSAFTModel, V, T, z,_data)
    return @f(a_hc,_data) + @f(a_disp,_data) + @f(a_mp,_data)
end

function a_res_minus_assoc(model::SAFTgammaMieModel, V, T, z,_data)
    dgc,X,vrdata = _data
    _,ρS,ζi,_ζ_X,_ζst,σ3x,m̄ = vrdata
    vrdata_disp = (dgc,ρS,ζi,_ζ_X,_ζst,σ3x,m̄)
    return @f(a_hs,_data) + a_disp(model,V,T,X,vrdata_disp)/sum(z) + @f(a_chain,_data)
end

function a_res_minus_assoc(model::SAFTVRMieModel, V, T, z, _data)
    return @f(a_hs,_data) + @f(a_dispchain,_data)
end

function a_res_minus_assoc(model ::SAFTVRMieGVModel, V, T, z, _data) 
    return @f(a_hs,_data) + @f(a_dispchain,_data) + @f(a_mp,_data)
end

function a_res_minus_assoc(model::SAFTVRSWModel, V, T, z, _data)
    return @f(a_mono,_data) + @f(a_chain,_data)
end

function a_res_minus_assoc(model::softSAFTModel, V, T, z, _data)
    return @f(a_LJ,_data) + @f(a_chain)
end

function a_res_minus_assoc(model::CPAModel,V,T,z,_data)
    return a_res(model.cubicmodel,V,T,z,_data)
end

#special dispatch for CPA
_X_and_Δ(model,V,T,z,neutral_data) = X_and_Δ(neutralmodel,V,T,z,neutral_data)

function _X_and_Δ(model::CPAModel,V0,T,z,neutral_data)
    n,ā,b̄,c̄ = _data
    V = V0 + c̄*n #volume translation
    return X_and_Δ(neutralmodel,V,T,z,neutral_data)
end

function a_res(model::ESElectrolyteModel,V,T,z,dep::DependentIonModel{MM1}) 
    rsp = dep.model
    neutralmodel = model.neutralmodel
    ionmodel = model.ionmodel
    neutral_data = data(neutralmodel,V0,T,z)
    Z = model.charge
    X,Δ = _X_and_Δ(neutralmodel,V,T,z,neutral_data)
    a_neutral = a_res_minus_assoc(neutralmodel,V,T,z,neutral_data) + a_assoc_impl(model,V,T,z,X)
    ϵ_r = __dielectric_constant(model, V, T, z, rsp, X, Δ)
    σ = get_sigma(ionmodel, V, T, z, neutralmodel, neutral_data)
    iondata = (Z, σ, ϵ_r)
    a_ion = a_res(ionmodel, V, T, z, iondata, neutralmodel, neutral_data)
    return a_neutral + a_ion
end

function dielectric_constant(model::ESElectrolyteModel, V, T, z,dep::DependentIonModel{MM1})
    neutralmodel = model.neutralmodel
    neutral_data = data(neutralmodel,V,T,z)
    X,Δ = X_and_Δ(neutralmodel,V,T,z,neutral_data)
    return __dielectric_constant(model, V, T, z, dep.model, X, Δ)
end

function __dielectric_constant(model::ESElectrolyteModel, V, T, z, RSPmodel::MM1, X, Δ)
    _1 = oneunit(Base.promote_eltype(RSPmodel,model.neutralmodel,V,T,z))
    _0 = zero(_1)
    assocmodel = model.neutralmodel
    μ0 = RSPmodel.params.mu.values
    γ = RSPmodel.params.gamma.values
    θ = RSPmodel.params.theta.values
    α = RSPmodel.params.polarizability.values
    z̄ = RSPmodel.params.coordz.values
    Z = model.charge
    n_neutral = count(iszero,Z)
    model_sites = getsites(assocmodel)
    sites = model_sites.i_sites
    ∑z = sum(z)
    ρ = N_A/V
    A = ρ/(3*ϵ_0)*sum(z[i]*α[i] for i ∈ @comps)
    ϵ_inf = (2*A+1)/(1-A)
    P = [fill(_0,n_neutral) for j ∈ 1:n_neutral]
    for i ∈ @ineutral
        if !isempty(sites[i])
            for j ∈ @ineutral
                if !isempty(sites[j])
                    #this should be a function
                    ∑ΔXX = sum(sum(Δ[i,j][a,b]*X[i][a]*X[j][b] for a ∈ sites[i]) for b ∈ sites[j])
                    P[i][j] = ρ/N_A*z[j]*∑ΔXX
                end
            end
        end
    end

    g = fill(_1,n_neutral)
    for i ∈ @ineutral
        μ0i = μ0[i]
        if !iszero(μ0[i])
            Pi = P[i]
            μ0i = μ0[i]
            ∑Pij = sum(Pi[j] for j ∈ @ineutral)
            gi = 1 + sum(z̄[i,j]*Pi[j]*cos(γ[i,j])/(∑Pij*cos(θ[i,j])+1)*μ0[j]/μ0i for j ∈ @ineutral)
            g[i] = gi
        end
    end

    B = ρ/(9*ϵ_0*k_B*T)*sum(z[i]*g[i]*μ0[i]^2 for i ∈ @ineutral)

    poly = (2,-(ϵ_inf+(ϵ_inf+2)^2*B),-ϵ_inf^2)

    return (-poly[2]+sqrt(poly[2]^2-8*poly[3]))/4
end
