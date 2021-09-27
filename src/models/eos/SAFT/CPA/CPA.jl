struct CPAParam <: EoSParam
    a::PairParam{Float64}
    b::PairParam{Float64}
    c1::SingleParam{Float64}
    Tc::SingleParam{Float64}
    epsilon_assoc::AssocParam{Float64}
    bondvol::AssocParam{Float64}
    Mw::SingleParam{Float64}
end

abstract type CPAModel <: SAFTModel end

struct CPA{T <: IdealModel,c <: CubicModel} <: CPAModel
    components::Array{String,1}
    icomponents::UnitRange{Int}
    cubicmodel::c
    params::CPAParam
    sites::SiteParam
    idealmodel::T
    absolutetolerance::Float64
    references::Array{String,1}
end

@registermodel CPA
export CPA

export CPA
function CPA(components; 
            idealmodel=BasicIdeal, 
            cubicmodel=RK, 
            alpha=CPAAlpha, 
            mixing=vdW1fRule,
            activity=nothing,
            translation=NoTranslation, 
            userlocations=String[], 
            ideal_userlocations=String[], 
            alpha_userlocations=String[],
            activity_userlocations=String[],
            mixing_userlocations=String[],
            translation_userlocations=String[],
            verbose=false)
    icomponents = 1:length(components)

    params,sites = getparams(components, ["SAFT/CPA", "properties/molarmass.csv","properties/critical.csv"]; userlocations=userlocations, verbose=verbose)
    Mw  = params["Mw"]
    k  = params["k"]
    Tc = params["Tc"]
    c1 = params["c1"]
    params["a"].values .*= 1E-1
    params["b"].values .*= 1E-3
    a  = epsilon_LorentzBerthelot(params["a"], k)
    b  = sigma_LorentzBerthelot(params["b"])
    epsilon_assoc = params["epsilon_assoc"]
    bondvol = params["bondvol"]
    packagedparams = CPAParam(a, b, c1, Tc, epsilon_assoc, bondvol,Mw)

    init_idealmodel = init_model(idealmodel,components,ideal_userlocations,verbose)
    init_alpha = init_model(alpha,components,alpha_userlocations,verbose)
    init_mixing = init_model(mixing,components,activity,mixing_userlocations,activity_userlocations,verbose)
    init_translation = init_model(translation,components,translation_userlocations,verbose)

    if occursin("RK",string(cubicmodel))
        cubicparams = RKParam(a, b, params["Tc"],params["pc"],Mw)
    elseif occursin("PR",string(cubicmodel))
        cubicparams = PRParam(a, b, params["Tc"],params["pc"],Mw)
    end

    init_cubicmodel = cubicmodel(components,icomponents,init_alpha,init_mixing,init_translation,cubicparams,init_idealmodel,1e-12,String[])

    references = ["10.1021/ie051305v"]

    model = CPA(components, icomponents, init_cubicmodel, packagedparams, sites, init_idealmodel, 1e-12, references)
    return model
end

function a_res(model::CPAModel, V, T, z)
    n = sum(z)
    ā,b̄,c̄ = cubic_ab(model.cubicmodel,V,T,z,n)
    return a_res(model.cubicmodel,V,T,z) + a_assoc(model,V+c̄*n,T,z)
end

ab_consts(model::CPAModel) = ab_consts(model.cubicmodel)

function a_assoc(model::CPAModel, V, T, z)
    x = z/∑(z)
    n_sites = model.sites.n_flattenedsites
    X_ = @f(X)
    return ∑(x[i]*∑(n_sites[i][a] * (log(X_[i][a])+(1-X_[i][a])/2) for a in @sites(i)) for i in @comps)
end

function X(model::CPAModel, V, T, z)
    _1 = one(V+T+first(z))
    x = z/∑(z)
    ρ = ∑(z)/V
    n_sites = model.sites.n_flattenedsites
    itermax = 100
    dampingfactor = 0.5
    error = 1.
    tol = model.absolutetolerance
    iter = 1
    X_ = [[_1 for a ∈ @sites(i)] for i ∈ @comps]
    X_old = deepcopy(X_)
    while error > tol
        iter > itermax && error("X has failed to converge after $itermax iterations")
        for i ∈ @comps, a ∈ @sites(i)
            rhs = 1/(1+∑(ρ*x[j]*∑(n_sites[j][b]*X_old[j][b]*@f(Δ,i,j,a,b) for b ∈ @sites(j)) for j ∈ @comps))
            X_[i][a] = (1-dampingfactor)*X_old[i][a] + dampingfactor*rhs
        end
        error = sqrt(∑(∑((X_[i][a] - X_old[i][a])^2 for a ∈ @sites(i)) for i ∈ @comps))
        X_old = deepcopy(X_)
        iter += 1
    end
    return X_
end

function Δ(model::CPAModel, V, T, z, i, j, a, b)
    ϵ_associjab = model.params.epsilon_assoc.values[i,j][a,b] * 1e2/R̄
    βijab = model.params.bondvol.values[i,j][a,b] * 1e-3
    x = z/∑(z)
    b = model.params.b.values
    b̄ = ∑(b .* (x * x'))
    η = b̄*∑(z)/(4*V)
    g = (1-η/2)/(1-η)^3
    bij = (b[i,i]+b[j,j])/2
    return g*(exp(ϵ_associjab/T)-1)*βijab*bij
end
