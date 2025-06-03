using Clapeyron, NLsolve, Plots

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

vvmax  = maximum(vv2)*0.7
vvmin  = minimum(vl2)*1.5

T1C = T1 .- 273.15
T2C = T2 .- 273.15
DH1 = (hV1.-hL1)./1e3
DH2 = (hV2.-hL2)./1e3

function A_critical(model,v,T,z)
    A(x) = Clapeyron.eos(model,x,T,z)
    dA(x) = Clapeyron.Solvers.derivative(A,x)
    d2A(x) = Clapeyron.Solvers.derivative(dA,x)
    d3A(x) = Clapeyron.Solvers.derivative(d2A,x)
    return -dA(v),-d2A(v),-d3A(v)
end

pc0, ∂p0_∂V, ∂²p0_∂V² = A_critical(model2,vc1,Tc1,[1.])

println("SAFTVRMieCP vc = $(vc1)")
println("SAFTVRMieCP pc = $(pc1)")
println("SAFTVRMieCP pc0 = $(pc0)")
println("SAFTVRMieCP ∂p0_∂V = $(∂p0_∂V)")
println("SAFTVRMieCP ∂²p0_∂V² = $(∂²p0_∂V²)")

A = zeros(3,3)
A[1,1] =  1.0/vc1^2
A[1,2] =  2.0/vc1^3
A[1,3] =  3.0/vc1^4
A[2,1] = -2.0/vc1^3
A[2,2] = -6.0/vc1^4
A[2,3] = -12.0/vc1^5
A[3,1] =  6.0/vc1^4
A[3,2] =  24.0/vc1^5
A[3,3] =  60.0/vc1^6
println("SAFTVRMieCP A = $(A)")
B = zeros(3)
B[1] =  pc1 - pc0
B[2] = -∂p0_∂V
B[3] = -∂²p0_∂V²
println("SAFTVRMieCP B = $(B)")
X = A \ B
println("SAFTVRMieCP X = $(X)")

pressure = zeros(N)
density = zeros(N)
for i in 1:N
    vv = vvmax + (vvmin - vvmax)*(i-1)/(N-1)
    density[i] = 1e-3/vv
    p,dpdV = Clapeyron.p∂p∂V(model2,vv,Tc1,[1.])
    pcrit = X[1]/vv^2 + 2.0*X[2]/vv^3 + 3.0*X[3]/vv^4
    pressure[i] = (p+pcrit)/1e5
end

p1 = plot([rhov1,rhol1,rhov2,rhol2,density], [psat1,psat1,psat2,psat2,pressure], xlabel = "rho [kmol/m3]", ylabel = "p [bara]",left_margin = 10Plots.mm,bottom_margin = 10Plots.mm,grid = :on,linewidth=3,size=(1600,1200))
display(p1)



