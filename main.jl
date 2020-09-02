# Initialize environment in current directory
#= @info("Ensuring example environment instantiated...") =#
#= import Pkg =#
#= Pkg.activate(@__DIR__) =#
#= Pkg.instantiate() =#

@info("Loading Zygote...")
using Zygote

include("constants.jl")

#### Some arbitrary conditions ####
temperature = 120
volume = 0.005
components = ["PMMA" "PS"]
compositions = [0.6 0.4]


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

model = PcSaft([1, 2], PcSaftParam(segments, sigmas, epsilons, ks)) # changed comp to [1, 2]
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
include("PcSaft.jl")
#= include("SaftGammaMie.jl") =#

#### Getting the gradiont ###
println(a_res(model, conditions))
function dif(conditions)
    global model
    return a_res(model, conditions)
end


gradients = gradient(dif, conditions)

println(gradients)


