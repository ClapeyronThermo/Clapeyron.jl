using Clapeyron

# a VLLE solution with a rich water phase and methanol partitioning between all phases

components = ["methane","ethane","propane","butane","hexane","octane","methanol","water"]
model = SAFTVRMie(components; assoc_options=AssocOptions(combining=:elliott))

p = 1.5e5
T = 15+273.15

zdry = [0.7,0.15,0.05,0.05,0.025,0.025]
xmeoh = 0.15
xwater = 0.05
z = append!(zdry*(1.0-xmeoh-xwater),xmeoh)
z = append!(z,xwater)

verbose = true
beta,xp,vp,Gsol = Clapeyron.tp_flash_impl(model,p,T,z, HELDTPFlash(verbose = verbose))

if !verbose
    for ip in eachindex(beta)
        println("Phase beta[$(ip)[] = $(beta[ip])")
    end
    println("Phase mole fraction:")
    for ip in eachindex(beta)
        println("Phase x[$(ip)[] = $(xp[ip])")
    end
    println("Phase volumes:")
    for ip in eachindex(beta)
        println("Phase volume[$(ip)[] = $(vp[ip])")
    end
    println("Minimum Gibbs Energy = $(Gsol)")
end