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
@newmodelsimple WalkerIdeal WalkerIdealModel WalkerIdealParam

export WalkerIdeal
function WalkerIdeal(components::Array{String,1}; userlocations::Array{String,1}=String[], verbose=false)
    params = getparams(components, ["ideal/WalkerIdeal.csv"]; userlocations=userlocations, verbose=verbose)
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
    return WalkerIdeal(packagedparams; references=references)
end

function a_ideal(model::WalkerIdealModel,V,T,z)
    x = z/sum(z)
    Mw = model.params.Mw.values
    Λ = @. h/√(k_B*T*Mw/N_A)
    Nrot = model.params.Nrot.values
    θ1 = model.params.theta1.values
    θ2 = model.params.theta2.values
    θ3 = model.params.theta3.values
    θ4 = model.params.theta4.values
    g1 = model.params.deg1.values
    g2 = model.params.deg2.values
    g3 = model.params.deg3.values
    g4 = model.params.deg4.values
    θ_vib = [θ1, θ2, θ2, θ2]
    g_vib = [g1, g2, g3, g4]
    return sum(x[i]*(log(z[i]*N_A/V*Λ[i]^3)-Nrot[i]/2*log(T)+sum(g_vib[v][i]*(θ_vib[v][i]/2/T+log(1-exp(-θ_vib[v][i]/T))) for v in 1:4)) for i in @comps)-1
end
