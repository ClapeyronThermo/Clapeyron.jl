#= using Plots =#
using JuliaSAFT

method = "PCSAFT"
model = system(["butane"], method)

p = 1e5
T = 500
z = create_z(model, [1.0])

a = JuliaSAFT.eos(model, z, p, T)
println(a)

#= (T_c, p_c, v_c) = get_Pcrit(model) =#
cp = get_isobaric_heat_capacity(model, z, p, T)
println(cp)


μ = get_chemical_potential(model, z, p, T)
println(μ)
g = get_Gibbs_free_energy(model, z, p, T)
println(g)
#= println(T_c) =#
#= temperature  = range(0.65*T_c,T_c,length=100) =#
#= (P_sat,v_l,v_v) = Psat(EoS,model,temperature) =#
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
