using Clapeyron, NLsolve, Plots

#components = ["carbon dioxide","nitrogen","water"]
components = ["carbon dioxide","water"]
model = GERG2008(components)

#p = 1.4e5*4.528*4.528
p = 15.25e5
T = 0.0+273.15

zdry = [1.0]
zwater=0.05
z=append!(zdry*(1-zwater),zwater)

pure = Clapeyron.split_pure_model(model)
crit = (Clapeyron.crit_pure).(pure)
vref = 0.0
for i = 1:length(z)
    Tc,pc,vc = crit[i]
    println("vc  = $(vc)")
    global vref += z[i]*vc
end

zs = [2.0176316286542775e-10, 1.0 - 2.0176316286542775e-10]

rho = Clapeyron.HELD_density(model,p,T,zs,vref)
mw = Clapeyron.molecular_weight(model,zs)
println("rho  = $(rho)")
vol = vref/rho
rhom = mw/vol
println("rhom  = $(rhom)")

#=
verbose = true
beta,xp,vp,Gsol = Clapeyron.tp_flash_impl(model,p,T,z, HELDTPFlash(verbose = verbose))

pure = Clapeyron.split_pure_model(model)
crit = (Clapeyron.crit_pure).(pure)
vref = 0.0
for i = 1:length(z)
    Tc,pc,vc = crit[i]
    global vref += z[i]*vc
end
=#
npoint = 1000
#comp = zeros(npoint)
density = zeros(npoint)
G = zeros(npoint)
dG = zeros(npoint)
ddG = zeros(npoint)
let rho = 0.15
    for i = 1:npoint
        global density[i] = rho
     #   z = [zCO2, 1.0 - zCO2]
     #   vol = Clapeyron.HELD_volume(model,p,T,z)
     #   rho = vref/vol
        Gi, dGi, ddGi = Clapeyron.HELD_volume2(model,p,T,zs,vref,rho)
        global G[i] = Gi
        global dG[i] = dGi
        global ddG[i] = ddGi
     #   zCO2 += (2 - 1)/npoint
        rho += (5.2 - 0.15)/npoint
    end
end

l = @layout [a b c]
p1 = plot(density, G, xlabel = "rho", ylabel = "G",left_margin = 10Plots.mm,bottom_margin = 10Plots.mm,grid = :on,linewidth=3,size=(1600,1200))
p2 = plot(density, dG, xlabel = "rho", ylabel = "dG",left_margin = 10Plots.mm,bottom_margin = 10Plots.mm,grid = :on,linewidth=3,size=(1600,1200))
p3 = plot(density, ddG, xlabel = "rho", ylabel = "ddG",left_margin = 10Plots.mm,bottom_margin = 10Plots.mm,grid = :on,linewidth=3,size=(1600,1200))
display(p1)
display(p2)
display(p3)
