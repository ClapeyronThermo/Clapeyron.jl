using Clapeyron, NLsolve

components = ["carbon dioxide","nitrogen","water"]
model = GERG2008(components)

p = 1.0e5
T = 20+273.15

zdry=[0.99, 0.01]
zwater=0.05
z=append!(zdry*(1-zwater),zwater)

verbose = false
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
println("Minimum Gibbs Energy = $(Gsol)")

zs = xp[1]
ps = p
pd = 4.0*ps
ts = T
td = 143.0+273.15
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

print("Compressor1 discharge temperature = ", round(tdis1.zero[1]-273.15,sigdigits=5)," deg C")

# add cooler to 30 deg C
tcooler1 = 30.0 + 273.15
betar,xpr,vpr,Gsol = Clapeyron.tp_flash_impl(model,pd,tcooler1,zs, HELDTPFlash(verbose = verbose))

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
println("Minimum Gibbs Energy = $(Gsol)")

zs = xp[1]
ps = pd
pd = 4.0*ps
ts = tcooler1
td = 156.0+273.15
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

print("Compressor2 discharge temperature = ", round(tdis2.zero[1]-273.15,sigdigits=5)," deg C")

# add cooler to 30 deg C
tcooler2 = 30.0 + 273.15
betar,xpr,vpr,Gsol = Clapeyron.tp_flash_impl(model,pd,tcooler2,zs, HELDTPFlash(verbose = verbose))

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
println("Minimum Gibbs Energy = $(Gsol)")

zs = xp[1]
ps = pd
pd = 4.0*ps
ts = tcooler2
td = 159.0+273.15
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

tdis3 = nlsolve(Compressor2 , [td])

print("Compressor3 discharge temperature = ", round(tdis3.zero[1]-273.15,sigdigits=5)," deg C")

# add cooler to 30 deg C
tcooler3 = 30.0 + 273.15
betar,xpr,vpr,Gsol = Clapeyron.tp_flash_impl(model,pd,tcooler3,zs, HELDTPFlash(verbose = verbose))

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
println("Minimum Gibbs Energy = $(Gsol)")
