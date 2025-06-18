using Clapeyron

# a VLLE solution with a rich water phase

components = ["methane","ethane","butane","hexane","octane","water"]
model = SAFTVRMie(components; assoc_options=AssocOptions(combining=:elliott))

p = 1.5e5
T = 15.0+273.15

zdry = [0.7,0.15,0.1,0.025,0.025]
xwater = 0.2
z = append!(zdry*(1.0-xwater),xwater)

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