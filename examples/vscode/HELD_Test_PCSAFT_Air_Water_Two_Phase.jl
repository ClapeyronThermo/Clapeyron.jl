using Clapeyron

components = ["nitrogen","oxygen","argon","carbon dioxide","water"]
model = PCSAFT(components; assoc_options=AssocOptions(combining=:elliott))

p = 1.5e5
T = 20+273.15

zdry=[0.7808,0.2095,0.0093,0.0004]
zwater=0.05
z=append!(zdry*(1-zwater),zwater)

isort = sortperm(z)
#is = [5,4,3,2,1]
println("Sorted model $(isort)")
models = Clapeyron.each_split_model(model, isort)
zs = Vector{Float64}(undef,length(z))
for i in eachindex(isort)
    zs[i] = z[isort[i]]
end
println("Sorted composition $(zs)")


verbose = true
betas,xps,vps,Gsol = Clapeyron.tp_flash_impl(models,p,T,zs, HELDTPFlash(verbose = verbose))

ivpsort = sortperm(vps, rev = true)
println("Sorted volumes $(ivpsort)")

println("Number phases found $(length(betas))")
beta = Vector{Float64}(undef,0)
for ip in eachindex(betas)
    push!(beta,betas[ivpsort[ip]])
end

xp = Vector{Vector{Float64}}(undef,0)
for ip in eachindex(betas)
    # undo sort
    xs = Vector{Float64}(undef,length(xps[ivpsort[ip]]))
    for i in eachindex(xps[ivpsort[ip]])
        xs[isort[i]] = xps[ivpsort[ip]][i]
    end 
    push!(xp,xs)
end

vp = Vector{Float64}(undef,0)
for ip in eachindex(vps)
    push!(vp,vps[ivpsort[ip]])
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
