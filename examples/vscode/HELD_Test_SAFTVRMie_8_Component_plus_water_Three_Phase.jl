using Clapeyron

components = ["methane","ethane","propane","butane","pentane","hexane","heptane","octane","water"]
model = SAFTVRMie(components; assoc_options=AssocOptions(combining=:elliott))

p = 5.0e5
T = 15.0+273.15

zdry = [0.7,0.15,0.025,0.025,0.025,0.025,0.025,0.025]
xwater = 0.1
z = append!(zdry*(1.0-xwater),xwater)

verbose = true

# note before you run this it takes hours to find the answer. However, its does succeed.
# speed is due to SAFT model and having associating components. Having nc = 9 is about as big as we can tolerate

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