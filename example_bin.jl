using Clapeyron

model = PCSAFT(["ethanol","water"])
toestimate = [
    Dict(
        :param => :epsilon,
        :indices => (1,2),
        :lower => model.params.epsilon.values[1,1]*0.9,
        :upper => model.params.epsilon.values[1,2]*1.1,
        :guess => model.params.epsilon.values[1,2]
    )
]

e = Estimation(model,toestimate,["bubble_comp.csv","bubble_temp.csv"])

# optimize!(e,Clapeyron.Metaheuristics.SA(x_initial=[e.toestimate.guess[1][1],e.toestimate.guess[2][1],e.toestimate.guess[3][1]],N=1000,tol_fun=1e-3))
optimize!(e,Clapeyron.Metaheuristics.ECA())

