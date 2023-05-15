#=
@setup_workload begin

    single = ["water"]
    pair = ["water", "ethanol"]
    @compile_workload begin
        s1 = PR(single)
        s2 = PCSAFT(single)
        s3 = SAFTVRMie(single)
        s4 = SAFTgammaMie(single)
        p1 = PCSAFT(pair)
        p2 = UNIFAC(pair)
        split_model(p1)
        split_model(p2)
        for si in (s1,s2,s3,s4)
            saturation_pressure(si,373.15)
            #crit_pure(si)
            pressure(si,0.03,373.15)
            volume(si,100000.0,373.15)
        end
        activity_coefficient(p2,100000,320,[0.5,0.5])
        #bubble_pressure(p1,320.0,[0.5,0.5])
        #bubble_pressure(p2,320.0,[0.5,0.5])
        #bubble_pressure(p1,320.0,Clapeyron.FractionVector(0.5))
        #bubble_pressure(p2,320.0,Clapeyron.FractionVector(0.5))
    end

end =#