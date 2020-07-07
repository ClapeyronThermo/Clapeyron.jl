#### Some arbitrary conditions ####
temperature = 298
volume = 1e-3
components = ["H2O" "CH4"]
compositions = [0.6 0.4]


#### Struct that contains all conditons ####
struct Conditions
    temperature::Float64
    volume::Float64
    components::Dict # contains components and compositons
end

conditions = Conditions(temperature, volume, Dict(zip(components, compositions)))


#### Struct that contains paramaters for PCSAFT ####
struct PcSaftParam
    segments::Dict
    sigmas::Dict
    epsilons::Dict
    ks::Dict # keys are tuples of every pair
end

struct SAFTVRMieParam
    segment::Dict
    sigma::Dict
    epsilon::Dict
    ShapeFactor::Dict
    lambdaA::Dict
    lambdaR::Dict
end


#### Types for SAFT models ####
abstract type Saft end

abstract type PcSaftFamily <: Saft end
abstract type SaftGammaMieFamily <: Saft end
abstract type SAFTVRMieFamily <: Saft end

struct PcSaft <: PcSaftFamily; components; parameters::PcSaftParam end
struct SPcSaft <: PcSaftFamily; components; parameters end
struct SaftGammaMie <: SaftGammaMieFamily; components; parameters end
struct SAFTVRMie <: SAFTVRMieFamily; components; parameters::SAFTVRMieParam end


#### Data from Pierre's script ####
import JSON

all_data = Dict()
open("all_data.json", "r") do f
    global all_data
    all_data = JSON.parse(f)  # parse and transform data
end
#### Some random test parameters ####
segment = Dict()
sigma = Dict()
epsilon = Dict()
ShapeFactor = Dict()
lambdaA = Dict()
lambdaR = Dict()

for k in components
    segment[k] = all_data["SEGMENT"][k]
    ShapeFactor[k] = all_data["SHAPEFACTOR"][k]

    #= for kk in intersect(keys(all_data["Binary_k"][k]), components) =#
    for kk in components
        sigma[(k,kk)] = all_data["SIGMA"][k][kk]*1e-10
        epsilon[(k,kk)] = all_data["EPSILON"][k][kk]
        lambdaA[(k,kk)] = all_data["LAMBDAA"][k][kk]
        lambdaR[(k,kk)] = all_data["LAMBDAR"][k][kk]
    end
end
model = SAFTVRMie(components, SAFTVRMieParam(segment, sigma, epsilon, ShapeFactor,lambdaA,lambdaR))

#### Functions relevant to PCSAFT ####
include("SAFTVRMie.jl")
#= include("SaftGammaMie.jl") =#


#### Calculation of a ####

println(a_mono(model, conditions))
