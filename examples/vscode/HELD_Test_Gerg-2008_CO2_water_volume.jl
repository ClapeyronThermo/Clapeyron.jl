using Clapeyron, Plots

components = ["carbon dioxide","nitrogen","water"]
model = GERG2008(components)
#model = SAFTVRMie(components)

#p = 64.0e5
#T = 30+273.15

#z=[0.9867980877492063, 0.009967685715362797, 0.0032342265354309285]
#zs=[0.00021287208725741366, 4.708156945730132e-9]
#zs=[0.4975080785711593, 0.0049838428576813995]
#z=append!(zs,1.0 - zs[1] - zs[2])

#vol = volume(model,p,T,z)
#mw = Clapeyron.molecular_weight(model,z)
#println("mw = $(mw)")
#rho = mw/vol
#println("Volume = $(vol) rho = $(rho)")
#rhom = mass_density(model,p,T,z)
#println("rhom = $(rhom)")

#verbose = true
#betar,xpr,vpr,Gsol = Clapeyron.tp_flash_impl(model,p,T,z, HELDTPFlash(verbose = verbose))
#=
density = zeros(250)
G = zeros(250)
dG = zeros(250)
let rho = 0.1
    for i = 1:250
        global density[i] = rho
        Gi, dGi = Clapeyron.HELD_volume2(model,p,T,z,rho)
        global G[i] = Gi
        global dG[i] = dGi
     #   println("rho = $(rho) G = $(Gi) and dG = $(dGi)")
        rho += 0.01
    end
end

l = @layout [a b]
p1 = plot(density, G, xlabel = "density", ylabel = "G",left_margin = 10Plots.mm,bottom_margin = 10Plots.mm,grid = :on,linewidth=3)
p2 = plot(density, dG, xlabel = "density", ylabel = "dG",left_margin = 10Plots.mm,bottom_margin = 10Plots.mm,grid = :on,linewidth=3)

plt = plot(p1, p2, layout = l, size=(1600,600))
display(plt)

print("Press any key to exit:")
n = readline()
=#
#G(x) = Clapeyron.HELD_volume(model,p,T,z,x)
#res = Clapeyron.Solvers.optimize(G,(0.1,1.0))
#res_min = Clapeyron.Solvers.x_minimum(res)
#println("rho = $(res_min)")
# Define objective function
# brent_f(x) = sin(x)
# Define the objective wrapper
#brent_scalar = Clapeyron.NLSolvers.ScalarObjective(; f = G)
# Define the optimization problem using the objective wrapper and bounds as a tuple
#brent_prob = Clapeyron.NLSolvers.OptimizationProblem(brent_scalar, (0.1, 1))
# Solve the problem using Brent's method for optimization
#Clapeyron.NLSolvers.solve(brent_prob, BrentMin(), OptimizationOptions())
#res = optimize(G, 0.001, 1.0, method = Brent())
#rho = Optim.minimizer(res)
#rho, found = Clapeyron.HELD_volume(model,p,T,z, 1e-10)
#println("bracket = $(rho) $(found)")

p = 1.0e5
T = 20+273.15

zdry=[0.99, 0.01]
zwater=0.05
z=append!(zdry*(1-zwater),zwater)

verbose = true
betar,xpr,vpr,Gsol = Clapeyron.tp_flash_impl(model,p,T,z, HELDTPFlash(verbose = verbose))

ivpsort = sortperm(vpr, rev = true)

println("Number phases found $(length(betar))")
beta = Vector{Float64}(undef,0)
for ip in eachindex(betar)
    push!(beta,betar[ivpsort[ip]])
end

xp = Vector{Vector{Float64}}(undef,0)
for ip in eachindex(betar)
    push!(xp,xpr[ivpsort[ip]])
end

vp = Vector{Float64}(undef,0)
for ip in eachindex(vpr)
    push!(vp,vpr[ivpsort[ip]])
end

for ip in eachindex(beta)
    println("Phase beta($(ip)) = $(beta[ip])")
end
println("Phase mole fraction:")
for ip in eachindex(beta)
    println("Phase x($(ip)) = $(xp[ip])")
end
println("Phase volumes:")
for ip in eachindex(beta)
    println("Phase volume($(ip)) = $(vp[ip])")
end
println("Phase densities:")
for ip in eachindex(beta)
    den = mass_density(model,p,T,xp[ip])
    println("Phase density($(ip)) = $(den)")
end
println("Minimum Gibbs Energy = $(Gsol)")

pure = Clapeyron.split_pure_model(model)
crit = (Clapeyron.crit_pure).(pure)
vref = 0.0
for i = 1:length(z)
    Tc,pc,vc = crit[i]
    global vref += z[i]*vc
end

println("vref = $(vref)")

for ip in eachindex(beta)
        println("Phase rho($(ip)) = $(vref/vp[ip])")
end

npoint = 100
density = zeros(npoint)
G = zeros(npoint)
dG = zeros(npoint)
ddG = zeros(npoint)
let rho = 1.6
    for i = 1:npoint
        global density[i] = rho
        Gi, dGi, ddGi = Clapeyron.HELD_volume2(model,p,T,xp[2],vref,rho)
        global G[i] = Gi
        global dG[i] = dGi
        global ddG[i] = ddGi
     #   println("rho = $(rho) G = $(Gi) and dG = $(dGi)")
        rho += 0.001
    end
end

l = @layout [a b c]
p1 = plot(density, G, xlabel = "density", ylabel = "G",left_margin = 10Plots.mm,bottom_margin = 10Plots.mm,grid = :on,linewidth=3,size=(1600,1200))
p2 = plot(density, dG, xlabel = "density", ylabel = "dG",left_margin = 10Plots.mm,bottom_margin = 10Plots.mm,grid = :on,linewidth=3,size=(1600,1200))
p3 = plot(density, ddG, xlabel = "density", ylabel = "ddG",left_margin = 10Plots.mm,bottom_margin = 10Plots.mm,grid = :on,linewidth=3,size=(1600,1200))
display(p1)
display(p2)
display(p3)
#plt = plot(p1, p2, p3, layout = l, size=(2400,600))
#display(plt)

Gr(x) = Clapeyron.HELD_func_rho(model,p,T,xp[2],vref,x)
dGr(x) = Clapeyron.Solvers.derivative(Gr,x)
ddGr(x) = Clapeyron.Solvers.derivative(dGr,x)

res = Clapeyron.Roots.find_zero(ddGr, 7.5)
println("Roots high side root = $(res)")

res = Clapeyron.Roots.find_zero(dGr, (res[1],7.5), no_pts = 15)
println("Roots = $(res)")

# calculate rho_ideal
rho_ideal = vref/(Clapeyron.R̄*T/p)
drho = rho_ideal/2.0

rho_min = 1.0e-6
rho_max = 7.5

if drho < rho_min
    rho_min = drho
    drho = 2.0*rho_min
else
    drho = (rho_min + rho_ideal)/2.0
end

rho1 = rho_min
rho2 = rho1+drho
if abs(dGr(rho2)) < sqrt(eps(Float64))
    rho2 += drho
end
rho_bracket = Vector{Vector{Float64}}(undef,0)
while rho2 <= rho_max
    if dGr(rho1)*dGr(rho2) < 0.0
        push!(rho_bracket,[rho1,rho2])
    end
    if rho1 >= 0.01
        global drho = 0.01
    end
    global rho1=rho2
    global rho2=rho1+drho
    if abs(dGr(rho2)) < sqrt(eps(Float64))
        global rho2 += drho
    end
end

println("Roots = $(rho_bracket)")

rho_found = Vector{Float64}(undef,0)
for ib = eachindex(rho_bracket)
    ans = Clapeyron.Roots.find_zero(dGr, rho_bracket[ib])
    for ia = eachindex(ans)
        stab = ddGr(ans[ia])
        if stab > 0.0
            push!(rho_found,ans[ia])
            println("Roots = $(ans[ia]) Gibbs = $(Gr(ans[ia]))")
        end
    end
end

rho_stable = rho_found[1]
if length(rho_found) > 1
    if Gr(rho_found[end]) < Gr(rho_stable)
        rho_stable = rho_found[end]
    end
end

println("Stable density = $(rho_stable)")

rho_test = Clapeyron.HELD_density(model,p,T,xp[2],vref)

println("Test density = $(rho_test)")

#print("Press any key to exit:")
#n = readline()

# So the HELD_func has multiple solution and the one in the middle is highly negative need to be able to find the real volume roots
# Try bracket using ddG = zero from above for liquid and below of vapout. Inbetween solution will have -ve ddG so an be discarded.
# once we have the liquid and vapour bracket we can use brent to minimise or indeed newton as we are close
# if we bracket vapout and the upper side is max denity then we have only 1 vapour root and likewise if we get min density for
# lower bracket for liquid we have only a liquid phase or its super critical and a fluid. Either was minimisation or newton will work.
