using Clapeyron

components = ["methane","ethane","propane","butane","water"]
model = PCSAFT(components; assoc_options=AssocOptions(combining=:elliott))

p = 90e5
T = 20+273.15

z = [0.64,0.06,0.28,0.005,0.015]

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