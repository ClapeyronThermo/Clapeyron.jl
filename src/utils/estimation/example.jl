using Clapeyron, BlackBoxOptim

model = UNIFAC(["methanol","benzene"])

toestimate = [
    Dict(
        :param => :A,
        :indices => (1,2),
        :symmetric => true,
        :lower => 120.,
        :upper => 200.,
        :guess => 150.
    )
]

e = Estimation(model,toestimate,["saturation_pressure.csv","saturation_liquid_density.csv"])

# optimize!(e,Clapeyron.Metaheuristics.SA(N=1000,tol_fun=1e-3))