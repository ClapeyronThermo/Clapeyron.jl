
# @info("Loading Zygote...")
# using Zygote
@info("Loading Optim...")
using Plots
using CSV, DataFrames
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

#### Some random test parameters ####
segment = Dict()
sigma = Dict()
epsilon = Dict()
lambdaA = Dict()
lambdaR = Dict()
components = ["Carbon dioxide"]

method = "SAFTVRMie"
include("Database.jl")
for k_ in components
    _, dF = lookup(k_, method)
    println(dF)
    i = 1
    segment[i] = dF[1,:m]
    sigma[i,i] = dF[1,:sigma]*1e-10
    epsilon[i,i] = dF[1,:epsilon]
    lambdaA[i,i] = dF[1,:lambdaA]
    lambdaR[i,i] = dF[1,:lambdaR]
    #= for kk in intersect(keys(all_data["Binary_k"][k]), components) =#

end

model = SAFTVRMie([1], SAFTVRMieParam(segment, sigma, epsilon, lambdaA,lambdaR)) # changed comp to [1, 2]
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
plt = plot(1000 ./temperature,P_sat,yaxis=:log,colour=:black)
plt = plot!([1000 ./T_c],[p_c],yaxis=:log,seriestype = :scatter,colour=:black)
plt = plot!(1000 ./T_exp,P_exp,yaxis=:log,seriestype = :scatter,colour=:red)

# plt = plot(temperature,H_vap)
