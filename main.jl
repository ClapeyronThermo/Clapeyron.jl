# Initialize environment in current directory
#= @info("Ensuring example environment instantiated...") =#
#= import Pkg =#
#= Pkg.activate(@__DIR__) =#
#= Pkg.instantiate() =#

@info("Loading Zygote...")
# using Zygote
@info("Loading Optim...")
using Optim,LinearAlgebra,Plots

include("constants.jl")

#### Some arbitrary conditions ####
temperature = 250
volume = 0.005
pressure = 1e6
components = ["CO2"]
compositions = [1]


#### Struct that contains all conditons ####
struct Conditions
    temperature::Float64
    volume::Float64
    compositions # A tuple
end

struct ReducedConditions
    temperature::Float64
    volume::Float64
    compositions # A list
end

#= conditions = Conditions(temperature, volume, Dict(zip(components, compositions))) =#
conditions = Conditions(temperature, volume, compositions)
#= conditions = Dict("temperature" => temperature, "volume" => volume, "components" => Dict(zip(components, compositions))) =#

#= struct conditions =#
#=     temperature =#
#=     volume =#
#=     compositions =#
#= end =#


#### Struct that contains paramaters for PCSAFT ####
struct PcSaftParam
    segments::Dict
    sigmas::Dict
    epsilons::Dict
    ks::Dict # keys are tuples of every pair
end


#### Types for SAFT models ####
abstract type Saft end

abstract type PcSaftFamily <: Saft end
abstract type SaftGammaMieFamily <: Saft end

struct PcSaft <: PcSaftFamily; components; parameters::PcSaftParam end
struct SPcSaft <: PcSaftFamily; components; parameters end
struct SaftGammaMie <: SaftGammaMieFamily; components; parameters end


#### Data from Pierre's script ####
import JSON

all_data = Dict()
open("all_data.json", "r") do f
    global all_data
    all_data = JSON.parse(f)  # parse and transform data
end


#### Some random test parameters ####
segments = Dict()
sigmas = Dict()
epsilons = Dict()
ks = Dict()

#= for k in components =#
#=     segments[k] = all_data["SEGMENT"][k] =#
#=     sigmas[k] = all_data["SIGMA"][k] =#
#=     epsilons[k] = all_data["EPSILON"][k] =#
#=     #1= for kk in intersect(keys(all_data["Binary_k"][k]), components) =1# =#
#=     for kk in components =#
#=             ks[(k,kk)] = k != kk ? all_data["Binary_k"][k][kk] : 0 =#
#=     end =#
#= end =#

for k_ in components
    if k_ == components[1]
        k = 1
    elseif k_ == components[2]
        k = 2
    end
    segments[k] = all_data["SEGMENT"][k_]
    sigmas[k] = all_data["SIGMA"][k_]
    epsilons[k] = all_data["EPSILON"][k_]
    #= for kk in intersect(keys(all_data["Binary_k"][k]), components) =#
    for kk_ in components
        if kk_ == components[1]
            kk = 1
        elseif kk_ == components[2]
            kk = 2
        end
            ks[(k,kk)] = k_ != kk_ ? all_data["Binary_k"][k_][kk_] : 0
    end
end

model = PcSaft([1], PcSaftParam(segments, sigmas, epsilons, ks)) # changed comp to [1, 2]
#= model = PcSaft(components, PcSaftParam(segments, sigmas, epsilons, ks)) =#

include("CombiningRules.jl")

function normal_to_reduced(model::PcSaftFamily, conditions)
    σ = combining_σ(model.parameters.sigmas)
    ϵ = combining_ϵ(model.parameters.epsilons)
    T = conditions.temperature
    v = conditions.volume
    x = conditions.compositions
    return ReducedConditions(T/ϵ, v/σ^3, x)
end

reduced_conditions = normal_to_reduced(model, conditions)

#### Functions relevant to PCSAFT ####
include("Ideal.jl")
include("PcSaft.jl")
include("Methods.jl")
#= include("SaftGammaMie.jl") =#

#### Getting the gradiont ###
EoS(z,v,T)    = a_ideal(model,z,v,T)+a_res(model,z,v,T)
# f(x) = (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2
# x0 = [0.0, 0.0]
# println(a_res(model, conditions))
#
# function g!(G,v)
#     global f
#     G[1] = gradient(f,v)[1]
# end
# function h!(H,v)
#     global f
#     H[1] = gradient(gradient(f,v))[1]
# end
# println(a_res(model,compositions,v0[1],temperature))
# println(gradient(f,v0)[1]*8.314*temperature)
pressure = range(1e6,5e6,length=100)
temperature  = ones(size(pressure))*250
z  = ones(size(pressure))
Vol = Volume(EoS,model,z,pressure,temperature)
plt = plot(Vol,pressure,xaxis=:log,yaxis=:log,fmt=:png)
P_exp=[1
1.1
1.2
1.3
1.4
1.5
1.6
1.7
1.785
1.785
1.8
1.9
2
2.1
2.2
2.3
2.4
2.5
2.6
2.7
2.8
2.9
3
3.1
3.2
3.3
3.4
3.5
3.6
3.7
3.8
3.9
4
4.1
4.2
4.3
4.4
4.5
4.6
4.7
4.8
4.9
5]
V_exp = [0.0018779
0.0016869
0.0015272
0.0013916
0.0012749
0.0011733
0.0010839
0.0010045
0.00094353
4.21E-05
4.21E-05
4.21E-05
4.20E-05
4.20E-05
4.20E-05
4.20E-05
4.20E-05
4.20E-05
4.19E-05
4.19E-05
4.19E-05
4.19E-05
4.19E-05
4.19E-05
4.18E-05
4.18E-05
4.18E-05
4.18E-05
4.18E-05
4.18E-05
4.17E-05
4.17E-05
4.17E-05
4.17E-05
4.17E-05
4.17E-05
4.17E-05
4.16E-05
4.16E-05
4.16E-05
4.16E-05
4.16E-05
4.16E-05]
plt = plot!(V_exp,P_exp*1e6,xaxis=:log,yaxis=:log,fmt=:png,seriestype = :scatter)
display(plt)
