struct EOS_LNG <: MultiFluidModel
    components::Vector{String}
    properties::MultiFluidPropertyParam
    ideal::MultiFluidIdealParam
    single::MultiFluidSingleParam
    pair::MultiFluidPairParam
    references::Vector{String}
end

@registermodel EOS_LNG

function EOS_LNG(components::Vector{String})
    I = getnames_gerg2008(components)

    params = getparams_eos_lng()
    Mw = each_split_model(params[:Mw],I)
    rhoc = each_split_model(params[:rhoc],I)
    vc = SingleParam("critical volume (dm3/mol)",rhoc.components,1 ./ rhoc.values)
    Tc = each_split_model(params[:Tc],I)
    pc = each_split_model(params[:pc],I)
    lb_v = each_split_model(params[:lb_v],I)
    properties = MultiFluidPropertyParam(Mw,rhoc,vc,Tc,pc,lb_v)
    n0 = each_split_model(params[:n0],I)
    theta = each_split_model(params[:theta],I)
    gamma_T = each_split_model(params[:gamma_T],I)
    gamma_v = each_split_model(params[:gamma_v],I)
    beta_T = each_split_model(params[:beta_T],I)
    beta_v = each_split_model(params[:beta_v],I)
    ideal = MultiFluidIdealParam(n0,theta,gamma_T,gamma_v,beta_T,beta_v)
    n = each_split_model(params[:n],I) |> pack_vectors
    t = each_split_model(params[:t],I) |> pack_vectors
    d = each_split_model(params[:d],I) |> pack_vectors
    c = each_split_model(params[:c],I) |> pack_vectors
    single = MultiFluidSingleParam(n,t,d,c)
    nij = each_split_model(params[:nij],I) |> pack_vectors
    tij = each_split_model(params[:tij],I) |> pack_vectors
    dij = each_split_model(params[:dij],I) |> pack_vectors
    Fij = each_split_model(params[:Fij],I)
    eta_ij = each_split_model(params[:eta_ij],I) |> pack_vectors
    epsilon_ij = each_split_model(params[:epsilon_ij],I) |> pack_vectors
    beta_ij = each_split_model(params[:beta_ij],I) |> pack_vectors
    gamma_ij = each_split_model(params[:gamma_ij],I) |> pack_vectors
    pair = MultiFluidPairParam(nij,tij,dij,Fij,beta_ij,gamma_ij,eta_ij,epsilon_ij)
    references = ["10.1021/je300655b"]
    return EOS_LNG(components,properties,ideal,single,pair,references)
end

export EOS_LNG
#nij tij dij ηij εij βij γij
#C1,C4


function getparams_eos_lng()
    params = getparams_gerg2008()
    γ_T =params[:gamma_T].values
    γ_v = params[:gamma_v].values
    β_T = params[:beta_T].values
    β_v = params[:beta_v].values
    lng_newpairs = [(1,6),(1,7),(1,8),(1,9)]
    _β_T = [0.9421, 0.9405, 0.9082, 0.886]
    _γ_T = [1.0307, 0.9917, 1.03884, 0.993]
    _β_v = [1.035, 1.0434, 1.02874, 1.023]
    _γ_v = [1.118, 1.143, 1.13209, 1.076]

    for (k,tup) in pairs(lng_newpairs)
        i,j = tup
        γ_T[i,j] =  _γ_T[k]
        γ_v[i,j] =  _γ_v[k]
        γ_T[j,i] =  _γ_T[k]
        γ_v[j,i] =  _γ_v[k]

        β_T[i,j] = _β_T[k] 
        β_v[i,j] = _β_v[k]
        β_T[j,i] = 1/_β_T[k] 
        β_v[j,i] = 1/_β_v[k]
    end

    nij = params[:nij].values
    tij = params[:tij].values
    dij = params[:dij].values
    eta_ij = params[:eta_ij].values
    epsilon_ij =params[:epsilon_ij].values
    beta_ij = params[:beta_ij].values
    gamma_ij = params[:gamma_ij].values
#C1,C4i
    nij[1,6] = [0.7588, -0.4386, -0.02273, 45.05, -2.291, -62.51, 33.32, -12.14]
    tij[1,6] = [1.02, 0.71, 1.57, 3.41, 2.12, 3.28, 3.37, 3.4]
    dij[1,6] = [1,2,3,1,1,1,2,1]
    eta_ij[1,6] = [1.34, 1.45, 0.96, 1.33, 1.9]
    epsilon_ij[1,6] = [0.59, 1.9, 0.87, 1.12, 1.43]
    beta_ij[1,6] = [1.07, 1.06, 1.11, 1.2, 1.23]
    gamma_ij[1,6] = [0.0, 0.0, 0.0, 0.0, 0.0]
#C1,C4i
    nij[1,7] = [0.9396, -0.1439, -0.1413, 35.32, -4.216, 59.17, -76.68, -41.39]
    tij[1,7] = [1.43, 0.3, 1.2, 3.1, 1.78, 3.36, 2.7, 3.7]
    dij[1,7] = [1,2,3,1,1,1,2,1]
    eta_ij[1,7] = [1.87, 1.05, 1.78, 1.19, 2.0]
    epsilon_ij[1,7] = [1.73, 0.78, 1.75, 1.84, 1.71]
    beta_ij[1,7] = [1.67, 1.76, 1.02, 1.76, 1.06]
    gamma_ij[1,7] = [0.0, 0.0, 0.0, 0.0, 0.0]

#C1C5
    nij[1,8] = [0.03711, -0.12154, 27.086, -13.614, -14.45, -0.46867]
    tij[1,8] = [1.54, 0.95, 0.47, 0.9, 0.22, 2.65]
    dij[1,8] = [2,3,1,1,1,2]
    eta_ij[1,8] = [0.6, 0.43, 0.7, 1.4]
    epsilon_ij[1,8] = [0.7, 0.61, 0.7, 0.46]
    beta_ij[1,8] = [0.916, 0.87, 0.86, 2.8]
    gamma_ij[1,8] = [0.5,0.5,0.5,0.5]
#C1C5i
    nij[1,9] =[0.051, -0.158, -67.49, -88.27, 154.9, 3.725]
    tij[1,9] = [0.2, 0.53, 1.79, 2.1, 2.0, 0.2]
    dij[1,9] = [2,3,1,1,1,2]
    eta_ij[1,9] = [0.64, 0.39, 0.48, 1.0]
    epsilon_ij[1,9] = [0.5, 0.5, 0.5, 0.5]
    beta_ij[1,9] = [1.56, 1.33, 1.46, 2.7]
    gamma_ij[1,9] = [0.0, 0.0, 0.0, 0.0]

    return params
end

"""
    EOS_LNG <: MultiFluidModel
    EOS_LNG(components::Vector{String})

## Imput Parameters

None

## Description

EOS-LNG: A Fundamental Equation of State for the Calculation of Thermodynamic Properties of Liquefied Natural Gases. valid for 21 compounds (`Clapeyron.GERG2008_names`). the EoS has new binary-specific parameters for methane + n-butane, methane + isobutane, methane + n-pentane, and methane + isopentane.

It uses the same functional form as [`GERG2008`](@ref).

## References

1. Thol, M., Richter, M., May, E. F., Lemmon, E. W., & Span, R. (2019). EOS-LNG: A fundamental equation of state for the calculation of thermodynamic properties of liquefied natural gases. Journal of Physical and Chemical Reference Data, 48(3), 033102. [doi:10.1063/1.5093800](https://doi.org/10.1063/1.5093800)
2. Kunz, O., & Wagner, W. (2012). The GERG-2008 wide-range equation of state for natural gases and other mixtures: An expansion of GERG-2004. Journal of Chemical and Engineering Data, 57(11), 3032–3091. [doi:10.1021/je300655b](https://doi.org/10.1021/je300655b)
"""
EOS_LNG