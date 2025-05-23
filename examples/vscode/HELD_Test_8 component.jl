using Clapeyron, Plots

components = ["methane","ethane","propane","butane","pentane","hexane","heptane","octane","water"]
model = SAFTVRMie(components; assoc_options=AssocOptions(combining=:elliott))

p = 5.0e5
T = 15+273.15

zdry = [0.883,0.08,0.021,0.007,0.005,0.002,0.001,0.001]
xwater = 0.1
zs = append!(zdry*(1.0-xwater),xwater)

pure = Clapeyron.split_pure_model(model)
crit = Clapeyron.crit_pure.(pure)
vref = 0.0
for i = 1:length(zs)
    Tc,pc,vc = crit[i]
    global vref += zs[i]*vc
end

z= [0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.992, 0.001, 0.001]

Vᵢ = Clapeyron.R̄*T/p
pᵢ,dpdVᵢ = Clapeyron.p∂p∂V(model,Vᵢ,T,z)
println("pᵢ,dpdVᵢ GERG2008 Clapeyron Volume method = $(pᵢ) $(dpdVᵢ)")
vol = Clapeyron.volume(model,p,T,z;threaded=false)
mw = Clapeyron.molecular_weight(model,z)
rhom = 1.0/vol*mw
println("Density Clapeyron Volume method = $(round(rhom,sigdigits=6)) kg/m3")

npoint = 1000
density = zeros(npoint)
#p1 = zeros(npoint)
#dp1 = zeros(npoint)
G1 = zeros(npoint)
dG1 = zeros(npoint)
ddG1 = zeros(npoint)
let rho = 0.2
    for i = 1:npoint
        global density[i] = rho/vref*mw 
        V = vref/rho
        #p1,dpdV1 = Clapeyron.p∂p∂V(model,V,T,z)
        #global dp1[i] = dpdV1
        Gi1, dGi1, ddGi1 = Clapeyron.HELD_volume2(model,p,T,z,vref,rho)
        global G1[i] = Gi1
        global dG1[i] = dGi1
        global ddG1[i] = ddGi1
        rho += (1.3 - 0.2)/npoint
    end
end

l = @layout [a b c d]
p1 = plot([density], [G1], xlabel = "rho", ylabel = "G",left_margin = 10Plots.mm,bottom_margin = 10Plots.mm,grid = :on,linewidth=3,size=(1600,1200))
p2 = plot([density], [dG1], xlabel = "rho", ylabel = "dG",left_margin = 10Plots.mm,bottom_margin = 10Plots.mm,grid = :on,linewidth=3,size=(1600,1200))
p3 = plot([density], [ddG1], xlabel = "rho", ylabel = "ddG",left_margin = 10Plots.mm,bottom_margin = 10Plots.mm,grid = :on,linewidth=3,size=(1600,1200))
#p4 = plot([density], [dp1], xlabel = "rho", ylabel = "dpdV",left_margin = 10Plots.mm,bottom_margin = 10Plots.mm,grid = :on,linewidth=3,size=(1600,1200))
display(p1)
display(p2)
display(p3)
#display(p4)

rho = Clapeyron.HELD_density(model,p,T,z,vref)
#mw = Clapeyron.molecular_weight(model,z)
println("rho  = $(rho)")
vol = vref/rho
rhom = mw/vol

println("Density = $(round(rhom,sigdigits=6)) kg/m3")