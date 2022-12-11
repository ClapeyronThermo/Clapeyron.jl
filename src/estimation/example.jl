#this file is meant to be loaded after Clapeyron itself is loaded.
using Clapeyron, Metaheuristics

function saturation_p(model::EoSModel,T)
    sat = saturation_pressure(model,T)
    return sat[1]
end

function saturation_rhol(model::EoSModel,T)
    sat = saturation_pressure(model,T)
    return 1/sat[2]
end

model = SAFTVRMie(["methane"])

toestimate = [
    Dict(
        :param => :epsilon,
        :lower => 130.,
        :upper => 300.,
        :guess => 250.
    ),
    Dict(
        :param => :sigma,
        :factor => 1e-10,
        :lower => 3.4,
        :upper => 4.0,
        :guess => 3.7
    )
    ,
    Dict(
        :param => :lambda_r,
        :lower => 10.,
        :upper => 16.,
        :guess => 12.
    )
]

e = Estimation(model,toestimate,["saturation_pressure.csv","saturation_liquid_density.csv"])

optimize!(e,Clapeyron.Metaheuristics.ECA())