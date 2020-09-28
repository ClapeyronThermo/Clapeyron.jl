#= using Plots =#
using JuliaSAFT,Plots

# model = system(["methanol"], "ogSAFT")
# p = 1.01e5
# T = 298.15
# z = create_z(model, [1.0])
# a = JuliaSAFT.eos(model, z, 1e-3, T)
# println(a)
# #
# cp = get_isobaric_heat_capacity(model, z, p, T,"liquid")
# println(cp)
# #
# μ = get_chemical_potential(model, z, p, T,"liquid")
# println(μ)
# g = get_Gibbs_free_energy(model, z, p, T,"liquid")
# println(g)
#
# (T_c, p_c, v_c) = get_Pcrit(model)
# println(T_c)
# Temp  = range(300,stop = T_c,length = 100)
# (P_sat,v_l,v_v) = get_Psat(model,Temp)
# println(P_sat)
# plt = plot(1e-3 ./v_l, Temp,color=:purple)
# plt = plot!(1e-3 ./v_v, Temp,color=:purple)
#
# println(sum(z[i]*μ[i] for i in model.components))
model_1 = system(["methanol"], "ogSAFT")
model_2 = system(["methanol"], "PCSAFT")
model_3 = system(["methanol"], "sPCSAFT")
model_4 = system(["methanol"], "SAFTVRMie")

(T_c, p_c, v_c) = get_Pcrit(model_1)
println(T_c)
Temp  = range(0.7*T_c,stop = T_c,length = 100)
(P_sat,v_l,v_v) = get_Psat(model_1,Temp)
println(P_sat)
plt = plot(1e-3 ./v_l, Temp,color=:purple)
plt = plot!(1e-3 ./v_v, Temp,color=:purple)

(T_c, p_c, v_c) = get_Pcrit(model_2)
println(T_c)
Temp  = range(0.7*T_c,stop = T_c,length = 100)
(P_sat,v_l,v_v) = get_Psat(model_2,Temp)
println(P_sat)
plt = plot!(1e-3 ./v_l, Temp,color=:red)
plt = plot!(1e-3 ./v_v, Temp,color=:red)
#
(T_c, p_c, v_c) = get_Pcrit(model_3)
println(T_c)
Temp  = range(0.7*T_c,stop = T_c,length = 100)
(P_sat,v_l,v_v) = get_Psat(model_3,Temp)
println(P_sat)
plt = plot!(1e-3 ./v_l, Temp,color=:orange)
plt = plot!(1e-3 ./v_v, Temp,color=:orange)
#
(T_c, p_c, v_c) = get_Pcrit(model_4)
println(T_c)
Temp  = range(0.7*T_c,stop = T_c,length = 100)
(P_sat,v_l,v_v) = get_Psat(model_4,Temp)
println(P_sat)
plt = plot!(1e-3 ./v_l, Temp,color=:blue)
plt = plot!(1e-3 ./v_v, Temp,color=:blue)
#= plt = plot(1 ./v_l,temperature,color=:red) =#
#= plt = plot!(1 ./v_v,temperature,color=:red) =#
#= plt = plot!(1 ./[v_c],[T_c],color=:red,seriestype = :scatter) =#
#= # plt = plot!(v_v_exp,T_exp,color=:red,seriestype = :scatter) =#
#= # plt = plot!(v_l_exp,T_exp,color=:red,seriestype = :scatter) =#
#= (model,EoS) = System(["n-butane"],"SAFTVRMie") =#
#= (T_c,p_c,v_c)=Pcrit(EoS,model) =#
#= println(T_c) =#
#= temperature  = range(0.65*T_c,T_c,length=100) =#
#= (P_sat,v_l,v_v) = Psat(EoS,model,temperature) =#
#= plt = plot!(1 ./v_l,temperature,color=:blue) =#
#= plt = plot!(1 ./v_v,temperature,color=:blue) =#
#= plt = plot!(1 ./[v_c],[T_c],color=:blue,seriestype = :scatter) =#
#= (model,EoS) = System(["n-butane"],"ogSAFT") =#
#= (T_c,p_c,v_c)=Pcrit(EoS,model) =#
#= println(T_c) =#
#= temperature  = range(0.65*T_c,T_c,length=100) =#
#= (P_sat,v_l,v_v) = Psat(EoS,model,temperature) =#
#= plt = plot!(1 ./v_l,temperature,color=:green) =#
#= plt = plot!(1 ./v_v,temperature,color=:green) =#
#= plt = plot!(1 ./[v_c],[T_c],color=:green,seriestype = :scatter) =#

#
# z = range(0.000001,stop=0.999999,length=100)
# h = []
# h_E = []
# h_1 = get_enthalpy(model,create_z(model,[0.999999,0.000001]),p,T,"liquid")
# h_2 = get_enthalpy(model,create_z(model,[0.000001,0.999999]),p,T,"liquid")
# for i in 1:length(z)
#     append!(h,get_enthalpy(model,create_z(model,[z[i],1-z[i]]),p,T,"liquid"))
#     append!(h_E,h[i]-z[i]*h_1-(1-z[i])*h_2)
# end
# print(h_E)
# plt = plot(z,h_E/1e3)
