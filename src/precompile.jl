using PrecompileTools
using Preferences
"""
    precompile_clapeyron!(val = true)

Activates or Deactivates the precompilation workload of the Clapeyron.,jl

!!! compat "Julia 1.9"
This function requires at least Julia 1.9
"""
function precompile_clapeyron! end
export precompile_clapeyron!

function precompile_clapeyron!(val = true)
    Preferences.set_preferences!(Clapeyron, "precompile_workload" => val; force=true)
    @info "Clapeyron's precompilation workload has been set to $(info_color(string(val))). this change will take effect on the next julia session."
end

@setup_workload begin

    single = ["water"]
    pair = ["water", "ethanol"]
    @compile_workload begin
        s1 = PR(single)
        s2 = PCSAFT(single)
        s3 = SAFTVRMie(single)
        s4 = SAFTgammaMie(single)
        s5 = RK(single)
        p1 = PCSAFT(pair)
        p2 = UNIFAC(pair)
        p3 = NRTL(pair)
        p4 = Wilson(pair)
        split_model(p1)
        split_model(p2)
        for si in (s1,s2,s3,s4,s5)
            saturation_pressure(si,373.15)
            #crit_pure(si)
            pressure(si,0.03,373.15)
            volume(si,100000.0,373.15)
        end
        activity_coefficient(p2,100000,320,[0.5,0.5])
        activity_coefficient(p3,100000,320,[0.5,0.5])
        activity_coefficient(p4,100000,320,[0.5,0.5])
        #bubble_pressure(p1,320.0,[0.5,0.5])
        #bubble_pressure(p2,320.0,[0.5,0.5])
        #bubble_pressure(p1,320.0,Clapeyron.FractionVector(0.5))
        #bubble_pressure(p2,320.0,Clapeyron.FractionVector(0.5))
    end
end