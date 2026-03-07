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
    
    @testset "mean ionic approach" begin
        ionmodel = SAFTVREMie(["water","ethanol"],["sodium","chloride"]; assoc_options=AssocOptions(combining=:elliott))
        model = MeanIonicApproach(ionmodel)
        p = 1e5
        T = 298.15
        m = [6.]
        zsolv = [0.5,0.5]
        z = molality_to_composition(ionmodel, salts, m, zsolv)
        w = Clapeyron.salt_compositions(model,z)
        z2 = Clapeyron.ion_compositions(model,w)
        @test z ≈ z2
        
        #mixed salt test
        ion_comps = ["water","1-but","K+","Cl-","Na+","SO4-2"]
        tsalts = ["KCl" => ["K+" => 1,"Cl-" => 1],"Na2SO4" => ["Na+" => 2, "SO4-2"=>1],"NaCl" => ["Na+" => 1,"Cl-" => 1]]
        tcharges = [0,0,1,-1,1,-2]
        salt_test = Clapeyron.explicit_salt_param(ion_comps,tsalts,tcharges)
        wsalt = [0.4,0.4,0.11,0.09,0.0]
        wion = [0.4,0.4,0.11,0.11,0.18,0.09]
        @test Clapeyron.ion_compositions(salt_test,wsalt) ≈ wion
        @test Clapeyron.salt_compositions(salt_test,wion) ≈ wsalt
        
        #rows are salts, cols are ions
        salt_mat = [1 0 0 0 0; 
                    0 1 0 0 0; 
                    0 0 1 0 0;
                    0 0 1 0 1; 
                    0 0 0 2 1; 
                    0 0 0 1 0]
        @test salt_test.salt_mat ≈ hcat(salt_mat,tcharges)
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
        
        #conversion between salt and ion base
        salt_model = MeanIonicApproach(model)
        res0_salt = Clapeyron.salt_compositions(salt_model,res0)
        res0_ion = Clapeyron.ion_compositions(salt_model,res0_salt)
        @test res0_ion.volumes ≈ res0.volumes
        @test res0_ion.fractions ≈ res0.fractions
        @test res0_ion.compositions[1] ≈ res0.compositions[1]
        @test res0_ion.compositions[2] ≈ res0.compositions[2]
        
        x_test =  [0.5522792450346276, 0.012878645399699698, 0.2174210525290201, 0.2174210570366527]
        @test x_test ≈ res0.compositions[1] rtol = 1e-6
        #@test x_test ≈ res1.compositions[1] rtol = 1e-6
        #@test x_test ≈ res2.compositions[1] rtol = 1e-6
    end

    @printline
end
