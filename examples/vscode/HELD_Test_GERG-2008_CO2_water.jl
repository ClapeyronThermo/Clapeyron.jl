using Clapeyron, NLsolve

components = ["carbon dioxide","nitrogen","water"]
model = EOS_CG(components)

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

# add cooler to 25 deg C and dehydration zwater needs to be low enough so a liquid water phase does not form in the JT calc
zdry=[zs[1], zs[2]]
zdry= zdry./sum(zdry)
zwater=0.0005
zs=append!(zdry*(1-zwater),zwater)
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

ps = pd
pd = ps/3.145
ts = tcooler3
td = -7.3+273.15
mw = Clapeyron.molecular_weight(model,zs)
hs=Clapeyron.VT_enthalpy(model,vp[1],ts,zs)/mw

function Valve1(F,x)
    beta,xp,vp,Gsol = Clapeyron.tp_flash_impl(model,pd,x[1],zs, HELDTPFlash(verbose = verbose))
    hd = 0.0
    for ip in eachindex(beta)
        hd += beta[ip]*Clapeyron.VT_enthalpy(model,vp[ip],x[1],xp[ip])
    end
    hd /= mw
    F[1] = hd - hs
end

tdis4 = nlsolve(Valve1 , [td])

Valve1_Tout= tdis4.zero[1]

println("Valve JT temperature = $(round(Valve1_Tout-273.15,sigdigits=5)) deg C")

beta,xp,vp,Gsol = Clapeyron.tp_flash_impl(model,pd,Valve1_Tout,zs, HELDTPFlash(verbose = verbose))

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