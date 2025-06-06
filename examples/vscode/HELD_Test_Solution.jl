using Clapeyron, Plots

components = ["methane","ethane","butane","hexane","octane","water"]
model = SAFTVRMie(components; assoc_options=AssocOptions(combining=:elliott))
#model = PCSAFT(components; assoc_options=AssocOptions(combining=:elliott))

p = 1.5e5
T = 15.0+273.15

zdry = [0.7,0.15,0.1,0.025,0.025]
xwater = 0.2
z = append!(zdry*(1.0-xwater),xwater)

pure = Clapeyron.split_pure_model(model)
crit = Clapeyron.crit_pure.(pure)
vref = 0.0
for i= 1:nc
    Tc,pc,vc = crit[i]
    global vref += z[i]*vc
end

xmin_SAFTVRMie = [[0.7123075723603035, 0.15242887318688825, 0.09941114344831163, 0.0190353971622987, 0.005079684746758808, 0.007795363091550251], 
        [8.192524850797215e-6, 2.1420937763026423e-5, 2.0090638215930553e-5, 3.337355836293869e-6, 1.007626601548929e-6, 6.913717598431754], 
        [0.005824659877551368, 0.008137650547407978, 0.0800630768619008, 0.21679358847466523, 0.6888426423112323, 0.8347834661570411]]

xmin_PCSAFT = [[0.006316560408822408, 0.009800276220052624, 0.08697232224811914, 0.21297590112307122, 0.6793353558745194, 0.870097120000675], 
        [0.7121006298597391, 0.15233700184106236, 0.09914706816128548, 0.01905596292411196, 0.005080307422121417, 0.008134048551555096], 
        [2.0390469365075103e-5, 1.6165576569853166e-5, 1.202826341931248e-5, 2.3599536294497255e-6, 5.51437942300973e-7, 6.63981431352832]]

x = xmin_SAFTVRMie[2]

nc = length(x)

zx = x[1:nc-1]
zn = 1.0 - sum(zx)
zi = append!(zx,zn)
rho = x[nc]
vol = vref/rho
mw = Clapeyron.molecular_weight(model,zi)
rhom = mw/vol

println("zᵢ $(zi)")
println("rho $(rho)")
println("vol $(vol)")
println("mw $(mw)")
println("density $(round(rhom,sigdigits=6)) kg/m3")

rho = Clapeyron.HELD_density(model,p,T,zi,vref)
mw = Clapeyron.molecular_weight(model,zi)
vol = vref/rho
rhom = mw/vol

println("rho $(rho)")
println("vol $(vol)")
println("mw $(mw)")
println("density $(round(rhom,sigdigits=6)) kg/m3")

#=
vref = 5.594807453383915e-5

z=[1.0]

Vᵢ = Clapeyron.R̄*T/p
pᵢ,dpdVᵢ = Clapeyron.p∂p∂V(model1,Vᵢ,T,z)
println("pᵢ,dpdVᵢ GERG2008 Clapeyron Volume method = $(pᵢ) $(dpdVᵢ)")
vol = Clapeyron.volume(model1,p,T,z;threaded=false)
mw = Clapeyron.molecular_weight(model1,z)
rhom = 1.0/vol*mw
println("Density GERG2008 Clapeyron Volume method = $(round(rhom,sigdigits=6)) kg/m3")

rho = Clapeyron.HELD_density(model2,p,T,z,vref)
mw = Clapeyron.molecular_weight(model2,z)
println("rho  = $(rho)")
vol = vref/rho
rhom = mw/vol

#println("Density IAPWS95 = $(round(rhom,sigdigits=6)) kg/m3")
println("Density SAFTVRMie = $(round(rhom,sigdigits=6)) kg/m3")

rho = Clapeyron.HELD_density(model1,p,T,z,vref)
mw = Clapeyron.molecular_weight(model1,z)
println("rho  = $(rho)")
vol = vref/rho
rhom = mw/vol

#println("Density GERG2008 = $(round(rhom,sigdigits=6)) kg/m3")
println("Density GERG2008 = $(round(rhom,sigdigits=6)) kg/m3")

npoint = 1000
#comp = zeros(npoint)
density = zeros(npoint)
p1 = zeros(npoint)
dp1 = zeros(npoint)
G1 = zeros(npoint)
dG1 = zeros(npoint)
ddG1 = zeros(npoint)
p2 = zeros(npoint)
dp2 = zeros(npoint)
G2 = zeros(npoint)
dG2 = zeros(npoint)
ddG2 = zeros(npoint)
let rho = 0.1
    for i = 1:npoint
        global density[i] = rho
     #   z = [zCO2, 1.0 - zCO2]
     #   vol = Clapeyron.HELD_volume(model,p,T,z)
     #   rho = vref/vol
        V = vref/rho
        p1,dpdV1 = Clapeyron.p∂p∂V(model1,V,T,z)
        global dp1[i] = dpdV1
        p2,dpdV2 = Clapeyron.p∂p∂V(model2,V,T,z)
        global dp2[i] = dpdV2
        Gi1, dGi1, ddGi1 = Clapeyron.HELD_volume2(model1,p,T,z,vref,rho)
        Gi2, dGi2, ddGi2 = Clapeyron.HELD_volume2(model2,p,T,z,vref,rho)
        global G1[i] = Gi1
        global dG1[i] = dGi1
        global ddG1[i] = ddGi1
        global G2[i] = Gi2
        global dG2[i] = dGi2
        global ddG2[i] = ddGi2
     #   zCO2 += (2 - 1)/npoint
        rho += (3.2 - 0.1)/npoint
    end
end

l = @layout [a b c d]
p1 = plot([density,density], [G1,G2], xlabel = "rho", ylabel = "G",left_margin = 10Plots.mm,bottom_margin = 10Plots.mm,grid = :on,linewidth=3,size=(1600,1200))
p2 = plot([density,density], [dG1,dG2], xlabel = "rho", ylabel = "dG",left_margin = 10Plots.mm,bottom_margin = 10Plots.mm,grid = :on,linewidth=3,size=(1600,1200))
p3 = plot([density,density], [ddG1,ddG2], xlabel = "rho", ylabel = "ddG",left_margin = 10Plots.mm,bottom_margin = 10Plots.mm,grid = :on,linewidth=3,size=(1600,1200))
p4 = plot([density,density], [dp1,dp2], xlabel = "rho", ylabel = "dpdV",left_margin = 10Plots.mm,bottom_margin = 10Plots.mm,grid = :on,linewidth=3,size=(1600,1200))
display(p1)
display(p2)
display(p3)
display(p4)
=#