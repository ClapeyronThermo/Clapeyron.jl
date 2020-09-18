
# @info("Loading Zygote...")
# using Zygote
@info("Loading Optim...")
using Plots
@info("Loading Zygote...")
#### Struct that contains paramaters for PCSAFT ####
struct ogSAFTParam
    segment::Dict
    sigma::Dict
    epsilon::Dict
end

#### Types for SAFT models ####
abstract type Saft end
abstract type ogSAFTFamily <: Saft end
struct ogSAFT <: ogSAFTFamily; components; parameters::ogSAFTParam end


#### Data from Pierre's script ####
import JSON

all_data = Dict()
open("all_data_ogSAFT.json", "r") do f
    global all_data
    all_data = JSON.parse(f)  # parse and transform data
end
#### Some random test parameters ####
segments = Dict()
sigmas = Dict()
epsilons = Dict()
components = ["C4H10"]
for k_ in components
    if k_ == components[1]
        k = 1
    elseif k_ == components[2]
        k = 2
    end
    segments[k] = all_data["SEGMENT"][k_]
    sigmas[k,k] = all_data["SIGMA"][k_]
    epsilons[k,k] = all_data["EPSILON"][k_]
    #= for kk in intersect(keys(all_data["Binary_k"][k]), components) =#
end

model = ogSAFT([1], ogSAFTParam(segments, sigmas, epsilons)) # changed comp to [1, 2]
#= model = PcSaft(components, PcSaftParam(segments, sigmas, epsilons, ks)) =#
#### Functions relevant to PCSAFT ####
include("Ideal.jl")
include("ogSAFT.jl")
include("Methods.jl")
#= include("SaftGammaMie.jl") =#

#### Getting the gradiont ###
EoS(z,v,T)    = 8.314*z[1]*T*(a_ideal(model,z,v,T)+a_res(model,z,v,T))
# println(a_seg(model,[1],1e-3,298))
# println(d(model,[1],1e-3,298,1))
(T_c,p_c,v_c) = Pcrit(EoS,model)
println(T_c)
temperature  = range(220,430,length=100)
(P_sat,v_v,v_l) = Psat(EoS,model,temperature)
T_exp = [220
245
270
295
320
345
370
395
420]
P_exp = [0.0078045
0.030885
0.091481
0.22034
0.45624
0.84406
1.4343
2.2868
3.4897]
plt = plot(1e4 ./temperature,P_sat/1e3,yaxis=:log,colour=:black)
plt = plot!([1e4 ./T_c],[p_c]/1e3,seriestype = :scatter,colour=:black)
plt = plot!(1e4 ./T_exp,P_exp*1e3,yaxis=:log,seriestype = :scatter,colour=:red)
