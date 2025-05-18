using Clapeyron, NLsolve, Plots

#components = ["carbon dioxide","nitrogen","water"]
components = ["carbon dioxide","water"]
model = GERG2008(components)

#p = 1.4e5*4.528*4.528*3.145
p = 1.4e5*4.528*4.528
#T = 25.0+273.15
T = -15+273.15

#zdry=[0.9750249927725483, 0.00185607409615914]
zdry = [1.0]
#zdry= zdry./sum(zdry)
zwater=0.05
z=append!(zdry*(1-zwater),zwater)
#z=[0.9750249927725483, 0.001856074096159142, 0.023118933131292534]
#z= z./sum(z)

verbose = true
beta,xp,vp,Gsol = Clapeyron.tp_flash_impl(model,p,T,z, HELDTPFlash(verbose = verbose))

hs = 0.0
for ip in eachindex(beta)
     global hs += beta[ip]*Clapeyron.VT_enthalpy(model,vp[ip],T,xp[ip])
end

println("hs = $(hs) J/mol")

pure = Clapeyron.split_pure_model(model)
crit = (Clapeyron.crit_pure).(pure)
vref = 0.0
for i = 1:length(z)
    Tc,pc,vc = crit[i]
    global vref += z[i]*vc
end

vol = Clapeyron.HELD_volume(model,p,T,z)
rho = vref/vol
println("rho = $(rho)")

npoint = 100
density = zeros(npoint)
G = zeros(npoint)
dG = zeros(npoint)
ddG = zeros(npoint)
let rho =1.7
    for i = 1:npoint
        global density[i] = rho
        Gi, dGi, ddGi = Clapeyron.HELD_volume2(model,p,T,z,vref,rho)
        global G[i] = Gi
        global dG[i] = dGi
        global ddG[i] = ddGi
     #   println("rho = $(rho) G = $(Gi) and dG = $(dGi)")
        rho += (1.72 - 1.7)/npoint
    end
end

l = @layout [a b c]
p1 = plot(density, G, xlabel = "density", ylabel = "G",left_margin = 10Plots.mm,bottom_margin = 10Plots.mm,grid = :on,linewidth=3,size=(1600,1200))
p2 = plot(density, dG, xlabel = "density", ylabel = "dG",left_margin = 10Plots.mm,bottom_margin = 10Plots.mm,grid = :on,linewidth=3,size=(1600,1200))
p3 = plot(density, ddG, xlabel = "density", ylabel = "ddG",left_margin = 10Plots.mm,bottom_margin = 10Plots.mm,grid = :on,linewidth=3,size=(1600,1200))
display(p1)
display(p2)
display(p3)