using Clapeyron

components = ["carbon dioxide","nitrogen","water"]
model = GERG2008(components)

p = 1.4e5*4.528*4.528*3.145
T = 25+273.15

zdry = [0.99,0.01]
zwater=0.0005
z=append!(zdry*(1-zwater),zwater)

verbose = false
beta,xp,vp,Gsol = Clapeyron.tp_flash_impl(model,p,T,z, HELDTPFlash(verbose = verbose))

println("Number phases found $(length(beta))")

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

zs = xp[1]
ps = p
pd = ps/3.145
ts = T
td = -8+273.15
mw = Clapeyron.molecular_weight(model,zs)
hs=Clapeyron.VT_enthalpy(model,vp[1],ts,zs)/mw

println("hs = $(hs)")

beta,xp,vp,Gsol = Clapeyron.tp_flash_impl(model,pd,td,zs, HELDTPFlash(verbose = verbose))

println("Number phases found $(length(beta))")

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

hd = 0.0
for ip in eachindex(beta)
    global hd += beta[ip]*Clapeyron.VT_enthalpy(model,vp[ip],td,xp[ip])
end
hd /= mw

println("hd = $(hd)")

println("hd - hs = $(hd - hs)")
