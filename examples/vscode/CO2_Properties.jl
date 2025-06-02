using Clapeyron, Plots

model1 = GERG2008(["carbon dioxide"])
model2 = SAFTVRMieCP(["carbon dioxide"])

(Tc1, pc1, vc1) = crit_pure(model1)
println("GERG2008 Tc, pc, rhoc = $(Tc1-273.15), $(pc1/1e5), $(1e-3/vc1)")
(Tc2, pc2, vc2) = crit_pure(model2)
println("SAFTVRMieCP Tc, pc, rhoc = $(Tc2-273.15), $(pc2/1e5), $(1e-3/vc2)")

N    = 100

T1    = LinRange(5+273.15, Tc1,  N)
T1C   = zeros(N)
psat1 = zeros(N)
vl1   = zeros(N)
vv1   = zeros(N)

hL1   = zeros(N)
hV1   = zeros(N)
cpL1  = zeros(N)
cpV1  = zeros(N)

v0 = [0.0,0.0]

for i in 1:N
    if i==1
        sat = saturation_pressure(model1, T1[i])
        psat1[i] = sat[1]/1e5
        vl1[i] = sat[2]
        vv1[i] = sat[3]
        global v0 = [vl1[i],vv1[i]]
    else
        sat = saturation_pressure(model1, T1[i]; v0=v0)
        psat1[i] = sat[1]/1e5
        vl1[i] = sat[2]
        vv1[i] = sat[3]
        global v0 = [vl1[i],vv1[i]]
    end
    hL1[i]  = Clapeyron.VT_enthalpy(model1,vl1[i],T1[i],[1.])
    hV1[i]  = Clapeyron.VT_enthalpy(model1,vv1[i],T1[i],[1.])
    cpL1[i] = Clapeyron.VT_isobaric_heat_capacity(model1,vl1[i],T1[i],[1.])
    cpV1[i] = Clapeyron.VT_isobaric_heat_capacity(model1,vv1[i],T1[i],[1.])
end

rhol1 = 1e-3 ./vl1
rhov1 = 1e-3 ./vv1

T2    = LinRange(5+273.15, Tc2,  N)
T2C   = zeros(N)
psat2 = zeros(N)
vl2   = zeros(N)
vv2   = zeros(N)

hL2   = zeros(N)
hV2   = zeros(N)
cpL2  = zeros(N)
cpV2  = zeros(N)

v0 = [0.0,0.0]

for i in 1:N
    if i==1
        sat = saturation_pressure(model2, T2[i])
        psat2[i] = sat[1]/1e5
        vl2[i] = sat[2]
        vv2[i] = sat[3]
        global v0 = [vl2[i],vv2[i]]
    else
        sat = saturation_pressure(model2, T2[i]; v0=v0)
        psat2[i] = sat[1]/1e5
        vl2[i] = sat[2]
        vv2[i] = sat[3]
        global v0 = [vl2[i],vv2[i]]
    end
    hL2[i]  = Clapeyron.VT_enthalpy(model2,vl2[i],T2[i],[1.])
    hV2[i]  = Clapeyron.VT_enthalpy(model2,vv2[i],T2[i],[1.])
    cpL2[i] = Clapeyron.VT_isobaric_heat_capacity(model2,vl2[i],T2[i],[1.])
    cpV2[i] = Clapeyron.VT_isobaric_heat_capacity(model2,vv2[i],T2[i],[1.])
end

rhol2 = 1e-3 ./vl2
rhov2 = 1e-3 ./vv2

T1C = T1 .- 273.15
T2C = T2 .- 273.15
DH1 = (hV1.-hL1)./1e3
DH2 = (hV2.-hL2)./1e3

p1 = plot([psat1,psat2], [T1C,T2C], xlabel = "p [bara]", ylabel = "T [oC]",left_margin = 10Plots.mm,bottom_margin = 10Plots.mm,grid = :on,linewidth=3,size=(1600,1200))
display(p1)
p2 = plot([rhov1,rhol1,rhov2,rhol2], [T1C,T1C,T2C,T2C], xlabel = "rho [kmol/m3]", ylabel = "T [oC]",left_margin = 10Plots.mm,bottom_margin = 10Plots.mm,grid = :on,linewidth=3,size=(1600,1200))
display(p2)
p3 = plot([T1C,T2C], [DH1,DH2], ylabel = "latent heat vaporisation [kJ/mol]", xlabel = "T [oC]",left_margin = 10Plots.mm,bottom_margin = 10Plots.mm,grid = :on,linewidth=3,size=(1600,1200))
display(p3)

