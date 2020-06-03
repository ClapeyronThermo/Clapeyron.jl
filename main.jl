#### Some arbitrary conditions ####
temperature = 120
volume = 0.005
components = ["PMMA" "PS"]
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

for k in components
    segments[k] = all_data["SEGMENT"][k]
    sigmas[k] = all_data["SIGMA"][k]
    epsilons[k] = all_data["EPSILON"][k]
    #= for kk in intersect(keys(all_data["Binary_k"][k]), components) =#
    for kk in components
            ks[(k,kk)] = k != kk ? all_data["Binary_k"][k][kk] : 0
    end
end
model = PcSaft(components, PcSaftParam(segments, sigmas, epsilons, ks))

#### Functions relevant to PCSAFT ####
include("PcSaft.jl")
#= include("SaftGammaMie.jl") =#


#### Calculation of a ####
function a_res(model::Saft, conditons)
    return a_hc(model, conditions) + a_disp(model, conditions)
end

println(a_res(model, conditions))
