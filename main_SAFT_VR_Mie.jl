
# @info("Loading Zygote...")
# using Zygote
@info("Loading Optim...")
# using Plots
@info("Loading Zygote...")
#### Struct that contains paramaters for PCSAFT ####
struct SAFTVRMieParam
    segment::Dict
    sigma::Dict
    epsilon::Dict
    lambdaA::Dict
    lambdaR::Dict
end

#### Types for SAFT models ####
abstract type Saft end
abstract type SAFTVRMieFamily <: Saft end
struct SAFTVRMie <: SAFTVRMieFamily; components; parameters::SAFTVRMieParam end


#### Data from Pierre's script ####
import JSON

all_data = Dict()
open("all_data_SaftVRMie.json", "r") do f
    global all_data
    all_data = JSON.parse(f)  # parse and transform data
end
#### Some random test parameters ####
segments = Dict()
sigmas = Dict()
epsilons = Dict()
lambdaAs = Dict()
lambdaRs = Dict()
components = ["CO2"]
for k_ in components
    if k_ == components[1]
        k = 1
    elseif k_ == components[2]
        k = 2
    end
    segments[k] = all_data["SEGMENT"][k_]
    sigmas[k,k] = all_data["SIGMA"][k_]
    epsilons[k,k] = all_data["EPSILON"][k_]
    lambdaAs[k,k] = all_data["LambdaA"][k_]
    lambdaRs[k,k] = all_data["LambdaR"][k_]
    #= for kk in intersect(keys(all_data["Binary_k"][k]), components) =#
end

model = SAFTVRMie([1], SAFTVRMieParam(segments, sigmas, epsilons, lambdaAs,lambdaRs)) # changed comp to [1, 2]
#= model = PcSaft(components, PcSaftParam(segments, sigmas, epsilons, ks)) =#
#### Functions relevant to PCSAFT ####
include("Ideal.jl")
include("SAFTVRMie.jl")
include("Methods.jl")
#= include("SaftGammaMie.jl") =#

#### Getting the gradiont ###
EoS(z,v,T)    = 8.314*z[1]*T*(a_ideal(model,z,v,T)+a_res(model,z,v,T))

temperature  = range(220,310,length=100)
z  = ones(size(temperature))
pressure = 101325*ones(size(temperature))
# (P_sat,v_v,v_l) = Psat(EoS,model,temperature)
H_vap=Enthalpy_of_vapourisation(EoS,model,temperature)
plt = plot(temperature,H_vap)
