#test structs for #154
abstract type PCSAFTModel_test <: SAFTModel end

# Defining the parameters used by the model
struct PCSAFTParam_test <: EoSParam
    Mw::SingleParam{Float64}
    segment::SingleParam{Float64}
    sigma::PairParam{Float64}
    epsilon::PairParam{Float64}
    epsilon_assoc::AssocParam{Float64}
    bondvol::AssocParam{Float64}
end

@testset "Reported errors" begin
    #https://github.com/ClapeyronThermo/Clapeyron.jl/issues/104
    @testset "#104" begin
        model = VTPR(["carbon dioxide"])
        p = 1e5
        T = 273.15
        @test fugacity_coefficient(model, p, T)[1] ≈ 0.9930275329424039 rtol = 1E-6
        @test activity_coefficient(model, p, T)[1] ≈ 1.0
    end

    @testset "#154" begin
        #there was a problem when using the @newmodel macros outside the Clapeyron module. this should suffice as a test.
        # Creating a model struct called PCSAFT, which is a sub-type of PCSAFTModel, and uses parameters defined in PCSAFTParam
        @newmodel PCSAFT_test PCSAFTModel_test PCSAFTParam_test
        @newmodelsimple PCSAFT_testsimple PCSAFTModel_test PCSAFTParam_test
        @newmodelgc PCSAFT_testgc PCSAFTModel_test PCSAFTParam_test

        #test if the macros actually generate something
        @test PCSAFT_test <: EoSModel #@newmodel
        @test PCSAFT_testsimple <: EoSModel #@newmodelsimple
        @test PCSAFT_testgc <: EoSModel #@newmodelgc
    end

    @testset "#162" begin
        #a longstanding problem, init_model didn't worked with functions.
        #a long time ago, SRK was a model, but now it is just a function that returns an RK model.
        model1 = Wilson(["water","ethanol"];puremodel=SRK)
        @test model1 isa Clapeyron.EoSModel

        #this case is just for compatibility with the notebooks that were originally released.
        model2 = VTPR(["carbon monoxide","carbon dioxide"];alpha = BMAlpha)
        @test model2 isa Clapeyron.EoSModel
    end

    @testset "##366" begin
        #=
        #366
        incorrect conversion of MixedGCSegmentParam.
        =#
        mix_segment_f64 = model.params.mixed_segment
        mix_segment_bigfloat = convert(Clapeyron.MixedGCSegmentParam{BigFloat},mix_segment_f64)
        @test mix_segment_f64.values.v ≈ mix_segment_bigfloat.values.v
    end

    @testset "#188" begin
        #=
        in cubics, we do a separation between the alpha model and the main model. while this allows for unprecedented
        flexibility, this also complicates the case of the default model. there is also an unrelated error, about needing to pass Vc,
        because our database has all the critical parameters in one file
        =#
        data = (
            species = ["A", "B"],
            Tc = [18.0, 3.5],
            Pc = [2.3, 0.1],
            Mw = [1.0, 1.0],
            acentricfactor = [0.1, 0.3]
        )

        file = ParamTable(
            :single,
            data,
            name="db1"
        )

        system = PR(["A", "B"], userlocations = [file])
        @test system.params.Tc[2] == 3.5
    end



    @testset "#212" begin
        #=
        Polar PCSAFT
        it uses `a` and `b` as site names
        =#
        model = PCPSAFT("water")
        @test model isa EoSModel
        @test "a" in model.sites.flattenedsites
        @test "b" in model.sites.flattenedsites
    end

    @testset "#75 - regression" begin
        #=
        pharmaPCSAFT, incorrect mixing rules
        =#

        userlocations = (Mw = [352.77,18.01],
        segment = [14.174,1.2046817736],
        sigma = [3.372,2.7927],
        epsilon = [221.261,353.9449],
        n_H = [2,1],
        n_e = [2,1],
        k = [0 -0.155;-0.155 0],
        kT = [0 2.641e-4;2.641e-4 0],
        epsilon_assoc = Dict(
        (("griseofulvin","e"),("griseofulvin","H")) => 1985.49,
        (("water08","e"),("water08","H")) => 2425.6714),
        bondvol = Dict(
        (("griseofulvin","e"),("griseofulvin","H")) => 0.02,
        (("water08","e"),("water08","H")) => 0.04509)
        )
        model = pharmaPCSAFT(["griseofulvin","water08"];userlocations = userlocations)
        #found while testing, test that we don't mutate the input
        @test userlocations.sigma == [3.372,2.7927]
        T = 310.15
        p = 1.01325e5
        z=[7.54e-7, 1-7.54e-7]

        γ1 = activity_coefficient(model,p,T,z)

        #this was an error too. check commit that added this
        #@test γ1[1] ≈ 55334.605821130834 rtol = 1e-4

        @test γ1[1] ≈ 51930.06908022231 rtol = 1e-4

    end


    @testset "SorptionModels.jl - init kij with user" begin
        #=
        on SL, passing k in userlocations did not work.
        =#


        v★(P★, T★,) = 8.31446261815324 * T★ / P★ / 1000000 # J / (mol*K) * K / mpa -> pa * m3 / (mol * mpa) ->  need to divide by 1000000 to get m3/mol
        ϵ★(T★) = 8.31446261815324 * T★ # J / (mol * K) * K -> J/mol
        r(P★, T★, ρ★, mw) = mw * (P★ * 1000000) / (8.31446261815324 * T★ * (ρ★ / 1e-6)) # g/mol * mpa * 1000000 pa/mpa / ((j/mol*K) * K * g/(cm3 / 1e-6 m3/cm3)) -> unitless

        P★ = [534., 630.]
        T★ = [755., 300.]
        ρ★ = [1.275, 1.515]
        mw = [100000, 44.01]
        kij = [0 -0.0005; -0.0005 0]
        model1 = Clapeyron.SL(
            ["PC", "CO2"],
            userlocations = Dict(
                :vol => v★.(P★, T★,),
                :segment => r.(P★, T★, ρ★, mw),
                :epsilon => ϵ★.(T★),
                :Mw => mw
            ),
            mixing_userlocations = (;k0 = kij, k1 = [0 0; 0 0], l = [0 0; 0 0])
        )

        model2 = Clapeyron.SL(
            ["PC", "CO2"],
            userlocations = Dict(
                :vol => v★.(P★, T★,),
                :segment => r.(P★, T★, ρ★, mw),
                :epsilon => ϵ★.(T★),
                :Mw => mw,
                :k => kij
            )
        )
        @test Clapeyron.get_k(model1)[1] ≈ Clapeyron.get_k(model2)[1]
    end

    @testset "https://github.com/ClapeyronThermo/Clapeyron.jl/discussions/239" begin
        #test for easier initialization of CPA/SAFT without association
        m3 = PCSAFT("water",userlocations =(segment = 1,Mw = 1,epsilon = 1,sigma = 1.0))
        @test length(m3.params.bondvol.values.values) == 0
    end

    @testset "infinite dilution derivatives - association" begin
        #reported by viviviena in slack
        #cp_params = (a = [36.54206320678348, 39.19437678197436, 25.7415, 34.91747774761048], b = [-0.03480434051958945, -0.05808483585041852, 0.2355, -0.014935581577635826], c = [0.000116818199785053, 0.0003501220208504329, 0.0001578, 0.000756101594841365], d = [-1.3003819534791665e-7, -3.6941157412454843e-7, -4.0939e-7, -1.0894144551347726e-6], e = [5.2547403746728466e-11, 1.276270011886522e-10, 2.1166e-10, 4.896983427747592e-10])
        #idealmodel = ReidIdeal(["water", "methanol", "propyleneglycol","methyloxirane"]; userlocations = cp_params)
        pcpsaft = PCPSAFT(["water", "methanol", "propyleneglycol","methyloxirane"])
        water = split_model(pcpsaft)[1]
        z = [1.0,0.0,0.0,0.0]
        v = volume(pcpsaft, 101325.0, 298.15, z, phase = :liquid)
        cp_pure = Clapeyron.VT_isobaric_heat_capacity(water, v, 298.15, 1.0)
        cp_mix = Clapeyron.VT_isobaric_heat_capacity(pcpsaft, v, 298.15, z)
        @test cp_mix ≈ cp_pure
        #@test cp_mix ≈ 69.21259493306137
    end

    @testset "#357 - electrolyte equilibria" begin

        model1 = ePCSAFT(["water"], ["calcium", "chloride"])
        salts1 = [("calcium chloride", ("calcium" => 1, "chloride" => 2))]
        x1 = molality_to_composition(model1, salts1, 1.0)
        bub_test = 374.7581484748338
        bubP1_test = 2971.744917038001

        bubT1 = bubble_temperature(model1, 101325, x1, FugBubbleTemperature(nonvolatiles = ["calcium", "chloride"]))[1]
        @test bubT1 ≈ bub_test rtol = 1e-6

        bubT1_chempot = bubble_temperature(model1, 101325, x1, ChemPotBubbleTemperature(T0=373.0, nonvolatiles=["calcium","chloride"]))[1]
        @test_broken bubT1_chempot ≈ bub_test rtol = 1e-6

        bubT1_chempot2 = bubble_temperature(model1, 101325, x1, ChemPotBubbleTemperature(nonvolatiles=["calcium","chloride"]))[1]
        @test_broken bubT1_chempot2 ≈ bub_test rtol = 1e-6

        bubT1_2 = bubble_temperature(model1, 101325, x1, FugBubbleTemperature(nonvolatiles = ["calcium", "chloride"],y0=[1.,0.,0.],vol0=(1.8e-5,1.),T0=373.15))[1]
        @test bubT1_2 ≈ bub_test rtol = 1e-6

        bubP1 = bubble_pressure(model1,298.15,x1,FugBubblePressure(nonvolatiles=["calcium","chloride"],y0=[1.,0.,0.],vol0=(1e-5,1.)))[1]
        @test bubP1 ≈ bubP1_test rtol = 1e-6

        bubP1_2 = bubble_pressure(model1,298.15,x1,FugBubblePressure(nonvolatiles=["calcium","chloride"]))[1]
        @test bubP1 ≈ bubP1_test rtol = 1e-6
    end

    @testset "#416" begin
        #if we build a cubic and only provide critical parameters (without acentric factor), read the database
        #to build the alpha model
        model = PR(["nitrogen"]; userlocations=(;
           Tc = [0.3],
           Pc = [1.],
           Mw = [1.]))
        @test !model.alpha.params.acentricfactor.ismissingvalues[1]
    end

    @testset "#513" begin
        #=
        error in x0_sat_pure_spinodal
        incorrect bounds, causing failure in convergence.
        =#
        model = tcRK("R1243zf")
        Tr = range(0.2,1.0,500)
        Tc = first(crit_pure(model))
        psat = first.(saturation_pressure.(model,Tr .* Tc))
        @test count(isnan,psat[100:end]) == 0
    end

    @testset "#528" begin
        #=
        Cubics: error in not considering if a and b are non-missing before building it 
        =#
        model = PR(["MethylLinoleate"];userlocations = (;a = [14.253080968127202], b = [0.0003800760318341279], Tc = [799.0001945392406], Pc = [1.3408208824579124e6], Mw = [294.472], Vc = [0.0012370135229711869] ))
        @test model.params.a.values[1,1] == 14.253080968127202
        @test model.params.b.values[1,1] == 0.0003800760318341279
    end

    @testset "#557" begin
        #=
        PCPSAFT: on the specific oxirane case, there were some points where the volume gets stuck between 2 iterations.
        solved via bracketing scheme in the volume solver.
        =#
        model = PCPSAFT(["water", "oxirane", "ethylene glycol"])
        TT = (85:0.1:100) .+ 273.15
        p = 17.5*101325
        v1 = volume.(model,p,TT,Ref([1.0,0.0,0.0]),phase = :l)
        v2 = volume.(model,p,TT,Ref([0.0,1.0,0.0]),phase = :l)
        v3 = volume.(model,p,TT,Ref([0.0,0.0,1.0]),phase = :l)
        @test 0 == count(isnan,v1)
        @test 0 == count(isnan,v2)
        @test 0 == count(isnan,v3)
    end
end
