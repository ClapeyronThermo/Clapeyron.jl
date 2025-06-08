using Clapeyron

components = ["nitrogen","oxygen","argon","carbon dioxide","water"]
model = PCSAFT(components; assoc_options=AssocOptions(combining=:elliott))

p = 1.5e5
T = 20+273.15

zdry=[0.7808,0.2095,0.0093,0.0004]
zwater=0.075
z=append!(zdry*(1-zwater),zwater)

verbose = true
beta,xp,vp,Gsol = Clapeyron.tp_flash_impl(model,p,T,z, HELDTPFlash(verbose = verbose))

if !verbose
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
end