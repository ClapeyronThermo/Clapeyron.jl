using Clapeyron, Plots

model = SAFTVRMie(["NO"];userlocations=(;
                    Mw = [30.0061],
                    epsilon=[167.93312110370104],
                    sigma=[2.6359936396800623],
                    shapefactor=[0.9364467128053224],
                    lambda_r=[45.75280817557093],
                    lambda_a=[6.0],
                    segment=[2.0],
                    n_H=[1.0],
                    n_e=[1.0],
                    epsilon_assoc = Dict((("NO","e"),("NO","H")) => 731.6569376461931),
                    bondvol = Dict((("NO","e"),("NO","H")) => 5.687465901480062e-29)))

v0 = Clapeyron.x0_crit_pure(model)
(Tc,pc,vc) = crit_pure(model,big.(v0))

N = 500

T = LinRange(Tc*0.7,Tc,N)

v0 = Clapeyron.x0_sat_pure(model,T[1])
p_sat = zeros(N)
rhol_sat = zeros(N)
rhov_sat = zeros(N)
h_vap = zeros(N)

for i in 1:N
    global v0, p_sat
    sat = saturation_pressure(model,T[i],v0=v0)
    v0 = [sat[2],sat[3]]
    p_sat[i] = sat[1]
    rhol_sat[i] = 1e-3/sat[2]
    rhov_sat[i] = 1e-3/sat[3]
    hl = Clapeyron.VT_enthalpy(model,sat[2],T[i])
    hv = Clapeyron.VT_enthalpy(model,sat[3],T[i])
    h_vap[i] = hv-hl
end

plt = plot(T,p_sat,yaxis=:log,label="")

savefig(plt,"pT.pdf")

plt = plot(rhol_sat,T,xaxis=:log,label="")
plot!(rhov_sat,T,xaxis=:log,label="")
savefig(plt,"rhoT.pdf")

plt = plot(T,h_vap,label="")

savefig(plt,"hvap.pdf")