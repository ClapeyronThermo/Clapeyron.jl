using Clapeyron, BlackBoxOptim

model = SAFTVRMie(["helium"])

toestimate = [
    Dict(
        :param => :epsilon,
        :lower => 3.7,
        :upper => 5.0,
        :guess => 3.8
    ),
    Dict(
        :param => :sigma,
        :factor => 1e-10,
        :lower => 3.3,
        :upper => 3.8,
        :guess => 3.5
    ),
    Dict(
        :param => :lambda_r,
        :lower => 12.0,
        :upper => 18.0,
        :guess => 16.0
    )
]

e = Estimation(model,toestimate,["saturation_pressure.csv","saturation_liquid_density.csv"])

optimize!(e)