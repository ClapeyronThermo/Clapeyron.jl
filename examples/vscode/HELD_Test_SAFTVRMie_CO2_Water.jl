using Clapeyron, NLsolve

components = ["carbon dioxide","nitrogen","water"]
model = SAFTVRMie(components; idealmodel=AlyLeeIdeal, assoc_options=AssocOptions(combining=:elliott))

p = 1.4e5
T = 25.0+273.15

zdry=[0.9981, 0.0019]
zwater=0.025
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
pd = 4.528*ps
ts = T
td = 161.0+273.15
mw = Clapeyron.molecular_weight(model,zs)
ds=mass_density(model,ps,ts,zs)
hs=enthalpy(model,ps,ts,zs)/mw
neta = 0.84

function Compressor1(F,x)
    ddis = mass_density(model,pd,x[1],zs)
    hdis = enthalpy(model,pd,x[1],zs)/mw
    npoly = log(pd/ps)/log(ddis/ds)
    F[1] = hdis - hs - npoly/(npoly-1)*(pd/ddis - ps/ds)/neta
end

tdis1 = nlsolve(Compressor1 , [td])

println("Compressor1 discharge temperature = $(round(tdis1.zero[1]-273.15,sigdigits=5)) deg C")

# add cooler to 25 deg C
tcooler1 = 25.0 + 273.15
beta,xp,vp,Gsol = Clapeyron.tp_flash_impl(model,pd,tcooler1,zs, HELDTPFlash(verbose = verbose))

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
ps = pd
pd = 4.528*ps
ts = tcooler1
td = 163.0+273.15
mw = Clapeyron.molecular_weight(model,zs)
ds=mass_density(model,ps,ts,zs)
hs=enthalpy(model,ps,ts,zs)/mw
neta = 0.84

function Compressor2(F,x)
    ddis = mass_density(model,pd,x[1],zs)
    hdis = enthalpy(model,pd,x[1],zs)/mw
    npoly = log(pd/ps)/log(ddis/ds)
    F[1] = hdis - hs - npoly/(npoly-1)*(pd/ddis - ps/ds)/neta
end

tdis2 = nlsolve(Compressor2 , [td])

println("Compressor2 discharge temperature = $(round(tdis2.zero[1]-273.15,sigdigits=5)) deg C")

# add cooler to 25 deg C
tcooler2 = 30.0 + 273.15
beta,xp,vp,Gsol = Clapeyron.tp_flash_impl(model,pd,tcooler2,zs, HELDTPFlash(verbose = verbose))

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
ps = pd
pd = 3.145*ps
ts = tcooler2
td = 136.0+273.15
mw = Clapeyron.molecular_weight(model,zs)
ds=mass_density(model,ps,ts,zs)
hs=enthalpy(model,ps,ts,zs)/mw
neta = 0.84

function Compressor3(F,x)
    ddis = mass_density(model,pd,x[1],zs)
    hdis = enthalpy(model,pd,x[1],zs)/mw
    npoly = log(pd/ps)/log(ddis/ds)
    F[1] = hdis - hs - npoly/(npoly-1)*(pd/ddis - ps/ds)/neta
end

tdis3 = nlsolve(Compressor3 , [td])

println("Compressor3 discharge temperature = $(round(tdis3.zero[1]-273.15,sigdigits=5)) deg C")

# add cooler to 30 deg C
tcooler3 = 25.0 + 273.15
beta,xp,vp,Gsol = Clapeyron.tp_flash_impl(model,pd,tcooler3,zs, HELDTPFlash(verbose = verbose))

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

components2 = ["carbon dioxide","nitrogen"]
model2 = SAFTVRMie(components2; idealmodel=AlyLeeIdeal, assoc_options=AssocOptions(combining=:elliott))

# remove water as temperature below freezing point of water
zs = [xp[1][1],xp[1][2]]
zs = zs./sum(zs)
#zs = xp[1]
ps = pd
pd = ps/3.145
ts = tcooler3
td = -7.4+273.15
mw = Clapeyron.molecular_weight(model2,zs)
hs=enthalpy(model2,ps,ts,zs)/mw

function Valve1(F,x)
    beta,xp,vp,Gsol = Clapeyron.tp_flash_impl(model2,pd,x[1],zs, HELDTPFlash(verbose = verbose))
    hdis = 0.0
    for ip in eachindex(beta)
    hdis += beta[ip]*enthalpy(model2,pd,x[1],xp[ip])
    end
    hdis /= mw
    F[1] = hdis - hs
end

tdis4 = nlsolve(Valve1 , [td])

Valve1_Tout= tdis4.zero[1]

print("Valve JT temperature = ", round(Valve1_Tout-273.15,sigdigits=5)," deg C")

beta,xp,vp,Gsol = Clapeyron.tp_flash_impl(model2,pd,Valve1_Tout,zs, HELDTPFlash(verbose = verbose))

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