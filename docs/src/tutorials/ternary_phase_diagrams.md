# Ternary phase diagrams
The notebook associated with this tutorial can be found [here](../../../examples/ternary_phase_diagrams.ipynb).

When considering three-component mixtures, it is now even more challenging to visualise large sections of the phase space on a 2D axis. Typically, ternary phase diagrams are visualised on a ternary plot, where the three axes correspond to the concentrations of the three components in equilibrium phases, at fixed temperature and pressure.

In order to actually plot these diagrams in Julia, one needs to install the ternary python package and call it from Julia:
```julia
julia> using PyCall

julia> ternary = pyimport("ternary")
```
## Liquid–Liquid Equilibrium
Ternary phase diagrams are typically used to visualise liquid–liquid equilibrium between two immiscible fluids and an entrainer. As an example, we'll consider a mixture of water, dichloromethane and acetone using UNIFAC:
```julia
julia> model = UNIFAC(["water",("dichloromethane",["CH2CL2"=>1]),("acetone",["CH3"=>1,"CH3CO"=>1])])
```
As we know water and dichloromethane will phase split, we can start our trace of the phase boundary on that side. To obtain the composition of the two phases, we need to use the `tp_flash(model, p, T, z)` method:
```julia
julia> p, T = 1e5, 298.15;

julia> z0 = [0.5,0.5,1e-10];

julia> (x,n,G) = tp_flash(model,p,T,z0,MichelsenTPFlash(equilibrium=:lle))
([0.9941813349915325 0.005818664997300019 1.1167552809179326e-11; 0.009492840315784318 0.9905071594960436 1.8817198921812175e-10], [0.49523586949546133 0.0028984768852708584 5.5629416193349775e-12; 0.004764130510680909 0.49710152310858685 9.443705837956089e-11], -6.933551486428005)
```
The difficulty with tracing the LLE region here is that we would ideally like to follow the Plait 
## Vapour–Liquid–Liquid Equilibrium
