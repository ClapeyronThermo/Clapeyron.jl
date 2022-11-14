using Clapeyron, BlackBoxOptim

function saturation_p(model::EoSModel,T)
    sat = saturation_pressure(model,T)
    return sat[1]
end

function saturation_rhol(model::EoSModel,T)
    sat = saturation_pressure(model,T)
    return 1/sat[2]
end

model = UNIFAC(["propane","benzene"])

toestimate = [
    Dict(
        :param => :A,
        :indices => (1,2),
        :lower => -1000.,
        :upper => 1000.,
        :guess => 0.
    )
]

e = Estimation(model,toestimate,["saturation_pressure.csv","saturation_liquid_density.csv"])

# optimize!(e,Clapeyron.Metaheuristics.ECA())