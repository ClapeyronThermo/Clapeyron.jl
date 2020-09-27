#= using Plots =#
using JuliaSAFT,Plots

method = "ogSAFT"
model = system(["methanol","butane"], "PCSAFT")

println(model)
p = 1e5
T = 300
z = create_z(model, [0.9,0.1])

a = JuliaSAFT.eos(model, z, 1e-3, T)
println(a)

cp = get_isobaric_heat_capacity(model, z, p, T,"liquid")
println(cp)
#
μ = get_chemical_potential(model, z, p, T,"liquid")
println(μ)
g = get_Gibbs_free_energy(model, z, p, T,"liquid")
println(g)

println(sum(z[i]*μ[i] for i in model.components))
#
# (T_c, p_c, v_c) = get_Pcrit(model_1)
# println(T_c)
# Temp  = range(0.7*T_c,stop = T_c,length = 100)
# (P_sat,v_l,v_v) = get_Psat(model_1,Temp)
# println(P_sat)
# plt = plot(1e-3 ./v_l, Temp,color=:purple)
# plt = plot!(1e-3 ./v_v, Temp,color=:purple)

# (T_c, p_c, v_c) = get_Pcrit(model_2)
# println(T_c)
# Temp  = range(0.7*T_c,stop = T_c,length = 100)
# (P_sat,v_l,v_v) = get_Psat(model_2,Temp)
# println(P_sat)
# plt = plot!(1e-3 ./v_l, Temp,color=:red)
# plt = plot!(1e-3 ./v_v, Temp,color=:red)
#
# (T_c, p_c, v_c) = get_Pcrit(model_3)
# println(T_c)
# Temp  = range(0.7*T_c,stop = T_c,length = 100)
# (P_sat,v_l,v_v) = get_Psat(model_3,Temp)
# println(P_sat)
# plt = plot!(1e-3 ./v_l, Temp,color=:orange)
# plt = plot!(1e-3 ./v_v, Temp,color=:orange)
#
# (T_c, p_c, v_c) = get_Pcrit(model_4)
# println(T_c)
# Temp  = range(0.7*T_c,stop = T_c,length = 100)
# (P_sat,v_l,v_v) = get_Psat(model_4,Temp)
# println(P_sat)
# plt = plot!(1e-3 ./v_l, Temp,color=:blue)
# plt = plot!(1e-3 ./v_v, Temp,color=:blue)
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
