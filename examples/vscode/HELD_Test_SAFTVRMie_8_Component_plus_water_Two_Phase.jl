using Clapeyron

components = ["methane","ethane","propane","butane","pentane","hexane","heptane","octane","water"]
model = SAFTVRMie(components; assoc_options=AssocOptions(combining=:elliott))

p = 5.0e5
T = 15.0+273.15

zdry = [0.883,0.08,0.021,0.007,0.005,0.002,0.001,0.001]
xwater = 0.1
z = append!(zdry*(1.0-xwater),xwater)

add_near_pure_guess = true
add_pure_component = fill(true,length(z))
add_random_guess = false
add_all_guess = false
verbose = true
Clapeyron.HELD_derivatives(model,p,T,z,add_near_pure_guess,add_pure_component,add_random_guess,add_all_guess,verbose)