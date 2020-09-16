
# @info("Loading Zygote...")
# using Zygote
@info("Loading Optim...")
using Plots
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

(T_c,p_c,v_c) = Pcrit(EoS,model)
temperature  = range(220,T_c,length=100)
(P_sat,v_v,v_l) = Psat(EoS,model,temperature)
# H_vap=Enthalpy_of_vapourisation(EoS,model,temperature)
T_exp = [220,230,240,250,260,270,280,290,300]
P_exp = [0.59913,0.89291,1.2825,1.785,2.4188,3.2033,4.1607,5.3177,6.7131]*1e6
plt = plot(1 ./temperature,P_sat,yaxis=:log)
plt = plot!([1 ./T_c],[p_c],yaxis=:log,seriestype = :scatter)
plt = plot!(1 ./T_exp,P_exp,yaxis=:log,seriestype = :scatter)

# plt = plot(temperature,H_vap)
