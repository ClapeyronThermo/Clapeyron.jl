using Clapeyron, Metaheuristics

method = ECA(;options=Options(iterations=100000))

model = SAFTVRMieCP(["carbon dioxide"];idealmodel=AlyLeeIdeal)

toestimate = [
    Dict(
        :param => :segment,
        :lower => 1.0,
        :upper => 2.0,
        :guess => 1.6936,
        :symmetric => false,
        :recombine => false
    ),
    Dict(
        :param => :dcrit,
        :lower => 10.0,
        :upper => 20.0,
        :guess => 15.5,
        :symmetric => false,
        :recombine => false
    ),
    Dict(
        :param => :epsilon,
        :lower => 100.,
        :upper => 300.,
        :guess => 207.89,
        :symmetric => true,
        :recombine => true
    ),
    Dict(
        :param => :sigma,
        :factor => 1e-10,
        :lower => 2.0,
        :upper => 4.0,
        :guess => 3.05,
        :symmetric => true,
        :recombine => true
    ),
    Dict(
        :param => :lambda_a,
        :lower => 22.0,
        :upper => 30.0,
        :guess => 26.408,
        :symmetric => true,
        :recombine => true
    ),
    Dict(
        :param => :lambda_r,
        :lower => 4.0,
        :upper => 6.0,
        :guess => 5.055,
        :symmetric => true,
        :recombine => true
    )
]

function sat_p_rv_rl_av_al(model::EoSModel,T)
    sat = saturation_pressure(model,T)
    p  = sat[1]
    vv = sat[3]
    vl = sat[2]
    rv = 1e-3/vv
    rl = 1e-3/vl
    av  = Clapeyron.VT_speed_of_sound(model,vv,T,[1.])
    al  = Clapeyron.VT_speed_of_sound(model,vl,T,[1.])
    return p,rv,rl,av,al
end

function sat_p_rl(model::EoSModel,T)
    sat = saturation_pressure(model,T)
    p  = sat[1]
    vl = sat[2]
    rl = 1e-3/vl
    return p,rl
end

estimator,objective,initial,upper,lower = Estimation(model,toestimate,["/home/cbranch/julia/dev/Clapeyron.jl/examples/data/sat_p_rl.csv"])

params, model = optimize(objective, estimator, method)

println("params = $(params)")