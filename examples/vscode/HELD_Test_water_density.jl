using Clapeyron, Plots

components = ["water"]
#model1 = EOS_CG(components)
model1 = GERG2008(components)
model2 = SAFTVRMie(components)
#model2 = IAPWS95()

p = 28e5
T = -7.5+273.15

vref = 5.594807453383915e-5

z=[1.0]

vol = Clapeyron.volume(model1,p,T,z)
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
G1 = zeros(npoint)
dG1 = zeros(npoint)
ddG1 = zeros(npoint)
G2 = zeros(npoint)
dG2 = zeros(npoint)
ddG2 = zeros(npoint)
let rho = 0.1
    for i = 1:npoint
        global density[i] = rho
     #   z = [zCO2, 1.0 - zCO2]
     #   vol = Clapeyron.HELD_volume(model,p,T,z)
     #   rho = vref/vol
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

l = @layout [a b c]
p1 = plot([density,density], [G1,G2], xlabel = "rho", ylabel = "G",left_margin = 10Plots.mm,bottom_margin = 10Plots.mm,grid = :on,linewidth=3,size=(1600,1200))
p2 = plot([density,density], [dG1,dG2], xlabel = "rho", ylabel = "dG",left_margin = 10Plots.mm,bottom_margin = 10Plots.mm,grid = :on,linewidth=3,size=(1600,1200))
p3 = plot([density,density], [ddG1,ddG2], xlabel = "rho", ylabel = "ddG",left_margin = 10Plots.mm,bottom_margin = 10Plots.mm,grid = :on,linewidth=3,size=(1600,1200))
display(p1)
display(p2)
display(p3)