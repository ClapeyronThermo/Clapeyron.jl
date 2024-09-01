using Clapeyron, PyCall
import PyPlot; const plt = PyPlot


# model1 = SAFTVRMieGV(["carbon dioxide"])
# model2 = GERG2008(["carbon dioxide"]) # using GERG as exp data
# models = [model1,model2]


# T_exp = range(190, 530,length=40)
# sat = saturation_pressure.(model2,T_exp)
# p_exp = [sat[i][1] for i ∈ 1:40] / 1e5
# v_l_exp = [sat[i][2] for i ∈ 1:40]
# v_v_exp = [sat[i][3] for i ∈ 1:40]

# T_mod = range(190, 530,length=100)
# sat = saturation_pressure.(model1,T_mod)
# p_mod = [sat[i][1] for i ∈ 1:100] / 1e5  # convert to bar
# v_l_mod = [sat[i][2] for i ∈ 1:100]
# v_v_mod = [sat[i][3] for i ∈ 1:100]


# plt.clf()
# plt.semilogx(1e-3 ./v_l_mod,T_mod,label="SAFT-VR Mie GV",linestyle=":",color="r")
# plt.semilogx(1e-3 ./v_v_mod,T_mod,label="",linestyle=":",color="r")
# plt.semilogx(1e-3 ./v_l_exp,T_exp,label="Experimental",marker="o",linestyle="",color="k")
# plt.semilogx(1e-3 ./v_v_exp,T_exp,label="",marker="o",linestyle="",color="k")
# plt.xlabel("Density / (mol/L)",fontsize=16)
# plt.ylabel("Temperature / K",fontsize=16)

# plt.semilogy(T_mod,p_mod,label="SAFT-VR Mie GV",linestyle=":",color="r")
# plt.semilogy(T_exp,p_exp,label="Experimental",marker="o",linestyle="",color="k")
# plt.xlabel("Temperature / K",fontsize=16)
# plt.ylabel("Pressure / bar",fontsize=16)
# plt.legend(loc="upper left",frameon=false,fontsize=12) 
# plt.xticks(fontsize=12)
# plt.yticks(fontsize=12)
# display(plt.gcf())


model = SAFTVRMieGV(["decane","heptane"])

p = 101325
T = 298.15
z = [0.03334, 0.96666]
Atest = helmholtz_free_energy_res(model,p,T,z)

# T = 313.15
# x = range(1e-5,1-1e-5,length=200)
# X = Clapeyron.FractionVector.(x)

# y = []
# p = []

# bub = bubble_pressure.(model,T,X)

# y = [bub[i][4][1] for i ∈ 1:200]
# p = [bub[i][1] for i ∈ 1:200]


# plt.clf()
# plt.plot(y,p./1e5,label="SAFT-VR Mie GV",linestyle="--",color="r")
# plt.plot(x,p./1e5,linestyle="--",color="r")

# plt.legend(frameon=false,fontsize=11) 
# plt.xlabel("composition of 1",fontsize=16)
# plt.ylabel("Pressure / bar",fontsize=16)
# plt.xticks(fontsize=12)
# plt.yticks(fontsize=12)
# plt.xlim([0,1])
# display(plt.gcf())