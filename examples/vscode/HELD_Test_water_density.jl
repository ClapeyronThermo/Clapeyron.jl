using Clapeyron, Plots

components = ["water"]
model = GERG2008(components)

p = 4.2e5
T = 20+273.15

z=[1.0]
mw = Clapeyron.molecular_weight(model,z)

vol = Clapeyron.HELD_volume(model,p,T,z)

rho = 1.0/vol*mw

println("Density = $(round(rho,sigdigits=5)) kg/m3")