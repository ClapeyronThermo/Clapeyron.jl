using Clapeyron

components = ["water","nitrogen","oxygen","argon","carbon dioxide"]
model = SAFTVRMie(components; assoc_options=AssocOptions(combining=:elliott))

p = 1.5e5
T = 20+273.15

#zdry=[0.7808,0.2095,0.0093,0.0004]
zdry=[0.0004,0.0093,0.2095,0.7808]
zwater=0.05
z=prepend!(zdry*(1-zwater),zwater)

verbose = true
add_all_guess = false
beta,xp,vp,Gsol = Clapeyron.tp_flash_impl(model,p,T,z, HELDTPFlash(add_all_guess = add_all_guess,verbose = verbose))

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