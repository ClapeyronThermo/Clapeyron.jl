#### Some arbitrary conditions ####
temperature = 120
volume = 0.005
components = ["PMMA" "PB"]
compositions = [0.6 0.4]


#### Struct that contains all conditons ####
struct Conditions
    temperature
    volume
    components::Dict # contains components and compositons
end

conditions = Conditions(temperature, volume, Dict(zip(components, compositions)))


#### Struct that contains paramaters for PCSAFT ####
struct PcSaftParam
    segments::Dict
    sigmas::Dict
    epsilons::Dict
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
    all_data=JSON.parse(f)  # parse and transform data
end


#### Some random test parameters ####
segments = Dict()
sigmas = Dict()
epsilons = Dict()

for k in components
    segments[k] = all_data["SEGMENT"][k]
    sigmas[k] = all_data["SIGMA"][k]
    epsilons[k] = all_data["EPSILON"][k]
end

model = PcSaft(components, PcSaftParam(segments, sigmas, epsilons))


#### Functions relevant to PCSAFT ####
include("PcSaft.jl")
#= include("SaftGammaMie.jl") =#


#### Calculation of a ####
function ares(model::Saft, conditons)
    return ahc(model, conditions) + adisp(model, conditions)
end

println(ares(model, conditions))
