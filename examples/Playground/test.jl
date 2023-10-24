using Clapeyron

model = SAFTgammaMie([("ammonia",["NH3"=>1])];userlocations=["params/."])

(Tc,pc,vc) = crit_pure(model)

(psat,vl,vv) = saturation_pressure(model,298.15)