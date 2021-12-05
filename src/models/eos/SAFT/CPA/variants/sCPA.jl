abstract type sCPAModel <: CPAModel end

struct sCPA{T <: IdealModel,c <: CubicModel} <: sCPAModel
    components::Array{String,1}
    icomponents::UnitRange{Int}
    cubicmodel::c
    params::CPAParam
    sites::SiteParam
    idealmodel::T
    assoc_options::AssocOptions
    absolutetolerance::Float64
    references::Array{String,1}
end

@registermodel sCPA
export sCPA

export sCPA
function sCPA(components; 
            idealmodel=BasicIdeal, 
            cubicmodel=RK, 
            alpha=sCPAAlpha, 
            mixing=vdW1fRule,
            activity=nothing,
            translation=NoTranslation, 
            userlocations=String[], 
            ideal_userlocations=String[], 
            alpha_userlocations=String[],
            activity_userlocations=String[],
            mixing_userlocations=String[],
            translation_userlocations=String[],
            verbose=false,
            assoc_options = AssocOptions())
    icomponents = 1:length(components)

    params,sites = getparams(components, ["SAFT/CPA/sCPA/", "properties/molarmass.csv","properties/critical.csv"]; userlocations=userlocations, verbose=verbose)
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

    model = sCPA(components, icomponents, init_cubicmodel, packagedparams, sites, init_idealmodel, assoc_options, 1e-12, references)
    return model
end

function Δ(model::sCPAModel, V, T, z, i, j, a, b)
    ϵ_associjab = model.params.epsilon_assoc.values[i,j][a,b] * 1e2/R̄
    βijab = model.params.bondvol.values[i,j][a,b] * 1e-3
    x = z/∑(z)
    b = model.params.b.values
    b̄ = ∑(b .* (x * x'))
    η = b̄*∑(z)/(4*V)
    g = (1-1.9η)^-1
    bij = (b[i,i]+b[j,j])/2
    return g*(exp(ϵ_associjab/T)-1)*βijab*bij/N_A
end