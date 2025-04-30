using Clapeyron

components = ["methane","ethane","propane","butane","hexane","octane","water"]
model = SAFTVRMie(components; assoc_options=AssocOptions(combining=:elliott))

p = 1.5e5
T = 15+273.15

zdry = [0.883,0.08,0.021,0.007,0.007,0.002]
xwater = 0.1
z = append!(zdry*(1.0-xwater),xwater)

verbose = true

beta,xp,vp,Gsol = Clapeyron.tp_flash_impl(model,p,T,z, HELDTPFlash(verbose = verbose))

if !verbose
    for ip = 1:length(beta)
        println("Phase beta($(ip)) = $(beta[ip])")
    end
    println("Phase mole fraction:")
    for ip = 1:length(beta)
        println("Phase x($(ip)) = $(xp[ip])")
    end
    println("Phase volumes:")
    for ip = 1:length(beta)
        println("Phase volume($(ip)) = $(vp[ip])")
    end
    println("Minimum Gibbs Energy = $(Gsol)")
end