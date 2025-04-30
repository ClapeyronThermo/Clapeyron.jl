using Clapeyron

components = ["methane","ethane","propane","butane","water"]
model = PCSAFT(components; assoc_options=AssocOptions(combining=:elliott))

p = 90e5
T = 20+273.15

z = [0.64,0.06,0.28,0.005,0.015]

add_near_pure_guess = true
add_pure_component = fill(true,length(z))
add_random_guess = false
add_all_guess = false
verbose = true
Clapeyron.HELD_derivatives(model,p,T,z,add_near_pure_guess,add_pure_component,add_random_guess,add_all_guess,verbose)