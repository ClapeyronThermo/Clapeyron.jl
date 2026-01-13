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

    system = ePCSAFT(["water","acetonitrile"],["sodium","chloride"])
    @testset "Tp flash" begin
        (x,n,G) = tp_flash(system,1e5,298.15,[0.4,0.4,0.1,0.1],MichelsenTPFlash(equilibrium=:lle,K0=[100.,1e-3,1000.,1000.]))
        @test_broken x[1,1] ≈ 0.07303405420993273 rtol = 1e-6
    end

    @printline
end
