struct WalkerIdealParam <: EoSParam
    Mw::SingleParam{Float64}
    Nrot::SingleParam{Int}
    theta1::SingleParam{Float64}
    theta2::SingleParam{Float64}
    theta3::SingleParam{Float64}
    theta4::SingleParam{Float64}
    deg1::SingleParam{Int}
    deg2::SingleParam{Int}
    deg3::SingleParam{Int}
    deg4::SingleParam{Int}
end

abstract type WalkerIdealModel <: IdealModel end
@newmodelgc WalkerIdeal WalkerIdealModel WalkerIdealParam

export WalkerIdeal
function WalkerIdeal(components::Array{String,1}; userlocations::Array{String,1}=String[], verbose=false)
    groups = GroupParam(components,["ideal/WalkerIdeal_Groups.csv"], verbose=verbose)
    params = getparams(groups, ["ideal/WalkerIdeal.csv"]; userlocations=userlocations, verbose=verbose)
    Mw = params["Mw"]
    Nrot = params["Nrot"]
    theta1 = params["theta_1"]
    theta2 = params["theta_2"]
    theta3 = params["theta_3"]
    theta4 = params["theta_4"]
    deg1 = params["deg_1"]
    deg2 = params["deg_2"]
    deg3 = params["deg_3"]
    deg4 = params["deg_4"]
    packagedparams = WalkerIdealParam(Mw, Nrot, theta1, theta2, theta3, theta4, deg1, deg2, deg3, deg4)
    references = ["10.1021/acs.jced.0c00723"]
    model = WalkerIdeal(packagedparams, groups, BasicIdeal, references=references, verbose=verbose)
    return model
end

function a_ideal(model::WalkerIdealModel,V,T,z)
    Mw = model.params.Mw.values
    Nrot = model.params.Nrot.values
    θ1 = model.params.theta1.values
    θ2 = model.params.theta2.values
    θ3 = model.params.theta3.values
    θ4 = model.params.theta4.values
    g1 = model.params.deg1.values
    g2 = model.params.deg2.values
    g3 = model.params.deg3.values
    g4 = model.params.deg4.values
    θ_vib = [θ1, θ2, θ3, θ4]
    g_vib = [g1, g2, g3, g4]
    n = model.groups.n_flattenedgroups
    res = zero(V+T+first(z))
    Σz = sum(z)
    x = z/Σz
    @inbounds for i in @comps
        ni = n[i]
        Mwi = sum(ni[k]*Mw[k] for k in @groups(i))
        Nroti = sum(ni[k]*Nrot[k] for k in @groups(i))/sum(ni[k] for k in @groups(i))
        Λ   = h/√(k_B*T*Mwi/N_A)
        res +=x[i]*(log(z[i]*N_A/V*Λ^3)-Nroti/2*log(T)+sum(ni[k]*sum(g_vib[v][k]*(θ_vib[v][k]/2/T+log(1-exp(-θ_vib[v][k]/T))) for v in 1:4) for k in @groups(i)))
    end
    return res - 1.
end
