struct CPAParam <: EoSParam
    epsilon_assoc::AssocParam{Float64}
    bondvol::AssocParam{Float64}
end

abstract type CPAModel <: EoSModel end

struct CPA{T <: IdealModel,c <: CubicModel} <: CPAModel
    components::Array{String,1}
    icomponents::UnitRange{Int}
    cubicmodel::c
    params::CPAParam
    sites::SiteParam
    idealmodel::T
    assoc_options::AssocOptions
    references::Array{String,1}
end
mw(model::CPAModel) = model.cubicmodel.params.Mw
molecular_weight(model::CPAModel,z=SA[1.0]) = comp_molecular_weight(mw(model),z)

@registermodel CPA
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
    verbose=false,
    assoc_options = AssocOptions())
    
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
    packagedparams = CPAParam(epsilon_assoc, bondvol)

    init_idealmodel = init_model(idealmodel,components,ideal_userlocations,verbose)
    init_alpha = init_model(alpha,components,alpha_userlocations,verbose)
    init_mixing = init_model(mixing,components,activity,mixing_userlocations,activity_userlocations,verbose)
    init_translation = init_model(translation,components,translation_userlocations,verbose)
    cubicparams = ABCubicParam(a, b, params["Tc"],params["pc"],Mw)

    init_cubicmodel = cubicmodel(components,icomponents,init_alpha,init_mixing,init_translation,cubicparams,init_idealmodel,String[])

    references = ["10.1021/ie051305v"]

    model = CPA(components, icomponents, init_cubicmodel, packagedparams, sites, init_idealmodel, assoc_options, references)
    return model
end


lb_volume(model::CPAModel,z = SA[1.0]) = lb_volume(model.cubicmodel,z)
T_scale(model::CPAModel,z=SA[1.0]) = T_scale(model.cubicmodel,z)
p_scale(model::CPAModel,z=SA[1.0]) = p_scale(model.cubicmodel,z)

function x0_crit_pure(model::CPAModel)
    lb_v = lb_volume(model)
    if isempty(model.params.epsilon_assoc.values[1,1])
        (2.0, log10(lb_v/0.3))
    else
        (2.75, log10(lb_v/0.3))
    end
end

data(model::CPAModel,V,T,z) = data(model.cubicmodel,V,T,z)

function a_res(model::CPAModel, V, T, z)
    _data = data(model,V,T,z)
    return a_res(model.cubicmodel,V,T,z,_data) + a_assoc(model,V,T,z,_data)
end

ab_consts(model::CPAModel) = ab_consts(model.cubicmodel)

function Δ(model::CPAModel, V, T, z, i, j, a, b,_data = data(model,V,T,z))
    Σz,ā,b̄,c̄ = _data
    ϵ_associjab = model.params.epsilon_assoc.values[i,j][a,b] * 1e2/R̄
    βijab = model.params.bondvol.values[i,j][a,b] * 1e-3
    V = V + c̄*Σz
    η = b̄/(4*V*Σz)
    g = (1-η/2)/(1-η)^3
    bij = (b[i,i]+b[j,j])/2
    return g*(exp(ϵ_associjab/T)-1)*βijab*bij/N_A
end