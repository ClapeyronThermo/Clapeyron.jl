using Pkg

Pkg.activate(".")
using Clapeyron, Measurements, ForwardDiffOverMeasurements
using PyCall
import PyPlot; const plt=PyPlot
corner = pyimport("corner")

N = 400

model = PCSAFT(["methane"];userlocations = (;
        Mw = [16.],
        epsilon = [150.03±9.39],
        sigma = [3.7039±0.27],
        segment = [1.0±0.08],
        k = [0.0;;],
        n_H=[0],
        n_e=[0],
        epsilon_assoc=nothing,
        bondvol=nothing),idealmodel=ReidIdeal)
system = GERG2008(["methane"])

(Tc1,pc1,vc1) = crit_pure(model)
(Tc2,pc2,vc2) = crit_pure(system)

T1 = LinRange(120.,Tc1.val,N)
T2 = LinRange(120.,Tc2,N)
T = LinRange(120.,300.,N)

sat_1 = saturation_pressure.(model,T1)
sat_2 = saturation_pressure.(system,T2)

cp = isobaric_heat_capacity.(model,4e6,T)
cp1 = [cp[i].val for i in 1:N]
cp1_err = [cp[i].err for i in 1:N]
cp2 = isobaric_heat_capacity.(system,4e6,T)

psat_1 = [sat_1[i][1].val for i in 1:N]
psat_err_1 = [sat_1[i][1].err for i in 1:N]
vl_1 = [sat_1[i][2].val for i in 1:N]
vl_err_1 = [sat_1[i][2].err for i in 1:N]
vv_1 = [sat_1[i][3].val for i in 1:N]
vv_err_1 = [sat_1[i][3].err for i in 1:N]


psat_2 = [sat_2[i][1] for i in 1:N]
vl_2 = [sat_2[i][2] for i in 1:N]
vv_2 = [sat_2[i][3] for i in 1:N]
vsat_2 = vcat(vl_2,reverse(vv_2))

# plt.clf()
# plt.semilogy(T2, psat_2,color="black")
# plt.plot(T1, psat_1,color="red")
# plt.plot(T1, psat_1,color="red")

# plt.fill_between(T1, psat_1-psat_err_1, psat_1+psat_err_1,alpha=0.2,color="red",edgecolor="white", linewidth=0.0)
# plt.ylabel(raw"$p$ / Pa")
# plt.xlabel(raw"$T$ / K")
# plt.xlim(minimum(T1),maximum(T1))
# plt.savefig("methane_psat_error.pdf")

# plt.clf()
# plt.plot(vl_1, T1,color="red")
# plt.plot(vv_1, T1,color="red")
# plt.fill_betweenx(T1,vl_1-vl_err_1, vl_1+vl_err_1,color="lightpink",edgecolor="white", linewidth=0.0)
# plt.fill_betweenx(T1,vv_1-vv_err_1, vv_1+vv_err_1,color="lightpink",edgecolor="white", linewidth=0.0)
# plt.semilogx(vsat_2, vcat(T2,reverse(T2)),color="black")
# plt.ylabel(raw"$T$ / K")
# plt.xlabel(raw"$v$ / (m$^3$/mol)")
# plt.ylim(minimum(T1),maximum(T1))
# plt.savefig("methane_rhosat_error.pdf")

# plt.clf()
# plt.plot(T, cp2,color="black")
# plt.plot(T, cp1,color="red")
# plt.fill_between(T,cp1-cp1_err, cp1+cp1_err,color="lightpink",edgecolor="white", linewidth=0.0)
# plt.ylim(0,300)
# plt.xlim(T[1],T[end])
# plt.ylabel(raw"$C_p$ / (J/K/mol)")
# plt.xlabel(raw"$T$ / K")

# plt.savefig("methane_cp_error.pdf")
