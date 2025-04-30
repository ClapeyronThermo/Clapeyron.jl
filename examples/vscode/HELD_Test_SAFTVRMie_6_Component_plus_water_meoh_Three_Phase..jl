using Clapeyron

components = ["methane","ethane","propane","butane","hexane","octane","methanol","water"]
model = SAFTVRMie(components; assoc_options=AssocOptions(combining=:elliott))

p = 1.5e5
T = 15+273.15

zdry = [0.7,0.15,0.05,0.05,0.025,0.025]
xmeoh = 0.15
xwater = 0.05
z = append!(zdry*(1.0-xmeoh-xwater),xmeoh)
z = append!(z,xwater)

add_near_pure_guess = true
add_pure_component = fill(true,length(z))
add_random_guess = false
add_all_guess = false
verbose = true
Clapeyron.HELD_derivatives(model,p,T,z,add_near_pure_guess,add_pure_component,add_random_guess,add_all_guess,verbose)