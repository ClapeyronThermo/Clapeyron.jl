using Clapeyron, Test

@testset "Electrolyte methods" begin
    @printline
    system = SAFTVREMie(["water"],["sodium","chloride"])
    salts = [("sodium chloride",["sodium"=>1,"chloride"=>1])]
    p = 1e5
    T = 298.15
    m = 1.
    z = molality_to_composition(system,salts,m)

    @testset "Bulk properties" begin
        @test mean_ionic_activity_coefficient(system,salts,p,T,m)[1] ≈ 0.6321355339999393 rtol = 1e-6
        @test osmotic_coefficient(system,salts,p,T,m)[1] ≈ 0.9301038212951828 rtol = 1e-6
        
        @test mean_ionic_activity_coefficient_sat(system,salts,T,m)[1] ≈ 0.632123986095978 rtol = 1e-6
        @test osmotic_coefficient_sat(system,salts,T,m)[1] ≈ 0.9300995995153606 rtol = 1e-5
    end

    fluid = SAFTVREMie(["water"],["sodium","chloride"])
    solid = SolidKs(["water","sodium.chloride","sodium.chloride.2water"])
    mapping = [(("water",1),)=>(("water",1)),
            (("sodium",1),("chloride",1))=>(("sodium.chloride",1)),
            (("sodium",1),("chloride",1),("water",2))=>(("sodium.chloride.2water",1))]

    system = CompositeModel(["water","sodium","chloride"];mapping=mapping,fluid=fluid,solid=solid)
    
    @testset "Solid solubility" begin
        @test sle_solubility(system,p,270.,[1,1,1];solute=["water"])[1] ≈ 0.9683019351679348 rtol = 1e-6
        @test sle_solubility(system,p,300.,[1,1,1];solute=["sodium.chloride"],x0=[-1.2])[1] ≈ 0.7570505178871523 rtol = 1e-6
    end

    @testset "electrolyte Tp flash" begin
        model = SAFTVREMie(["water","ethanol"],["sodium","chloride"]; assoc_options=AssocOptions(combining=:elliott))
        K0 = [1e3,1e-2,1000.,10.]
        method_1st_order = MichelsenTPFlash(equilibrium=:lle,K0=K0, nacc=0,ss_iters = 10)
        method_2nd_order = MichelsenTPFlash(equilibrium=:lle,K0=K0, nacc=0,ss_iters = 10,second_order = true)
        method_0th_order = RRTPFlash(equilibrium=:lle,K0=K0)
        
        p = 1e5
        T = 298.15
        m = [6.]
        zsolv = [0.5,0.5]
        z = molality_to_composition(model, salts, m, zsolv)

        res0 = Clapeyron.tp_flash2(model,p,T,z,method_0th_order)
        res1 = Clapeyron.tp_flash2(model,p,T,z,method_1st_order)
        res2 = Clapeyron.tp_flash2(model,p,T,z,method_2nd_order)

        x_test =  [0.5522792450346276, 0.012878645399699698, 0.2174210525290201, 0.2174210570366527]
        @test x_test ≈ res0.compositions[1] rtol = 1e-6
    end

    @printline
end
