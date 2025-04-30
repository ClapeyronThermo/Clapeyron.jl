using Clapeyron

components = ["hexane","aniline","1,1,2-trichloroethane","water","nitrogen"]
model = SRK(components;mixing=MHV2Rule, activity=UNIFAC)

p = 1.0e5
T = 5+273.15

z=ones(length(model))/length(model)

add_near_pure_guess = true
add_pure_component = fill(true,length(z))
add_random_guess = false
add_all_guess = false
verbose = true
Clapeyron.HELD_derivatives(model,p,T,z,add_near_pure_guess,add_pure_component,add_random_guess,add_all_guess,verbose)