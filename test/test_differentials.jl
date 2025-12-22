using Clapeyron, Test

@printline

struct TestModel <: EoSModel end
Clapeyron.idealmodel(::TestModel) = BasicIdeal()
function Clapeyron.eos_impl(model::TestModel,V,T,z)
    z_part = 1*z[1] + 2*z[2]+3*z[3]
    v_part = log(V)
    T_part = T^6
    vt_part = V*T*V*T
    return vt_part+T_part+v_part+z_part
end

function Clapeyron.∂f∂V(model::TestModel,V,T,z::AbstractVector)
    f(v) = eos(model,v,T,z)
    return Clapeyron.Solvers.derivative(f,V)
end

∂f∂T_analytical(model::TestModel,V,T,z) =6*T^5 +2*T*V*V
∂f∂V_analytical(model::TestModel,V,T,z) =1/V +2*V*T*T
∂2f∂T2_analytical(model::TestModel,V,T,z) =5*6*T^4 +2*V*V
∂2f∂V2_analytical(model::TestModel,V,T,z) =-1/(V*V) +2*T*T
∂2f∂V∂T_analytical(model::TestModel,V,T,z) =4*V*T
∂3f∂V3_analytical(model::TestModel,V,T,z) =2/(V*V*V)
∂3f∂V2∂T_analytical(model::TestModel,V,T,z) = 4*T
∂3f∂V∂T2_analytical(model::TestModel,V,T,z) = 4*V


@testset "differentials" begin
    model = TestModel()
    T = 500*rand()
    V = 10*rand()
    z = rand(3)
    f = Clapeyron.eos(model,V,T,z)
    v = ∂f∂V_analytical(model,V,T,z)
    t = ∂f∂T_analytical(model,V,T,z)
    vv = ∂2f∂V2_analytical(model,V,T,z)
    vt = ∂2f∂V∂T_analytical(model,V,T,z)
    tt = ∂2f∂T2_analytical(model,V,T,z)
    vvv = ∂3f∂V3_analytical(model,V,T,z)

    df = [v,t]
    d2f = [vv vt;vt tt]
    p = -v
    pv = -vv
    pvv = -vvv
    pt = -vt
    dp = [pv,pt]
    pvt = -∂3f∂V2∂T_analytical(model,V,T,z)
    ptt = -∂3f∂V∂T2_analytical(model,V,T,z)
    d2p = [pvv pvt;pvt ptt]

    @testset "first order" begin
        @test Clapeyron.∂f∂T(model,V,T,z) ≈ t
        @test Clapeyron.∂f∂V(model,V,T,z) ≈ v
        df1 = Clapeyron.∂f(model,V,T,z)
        @test df1[2] ≈ f
        @test df1[1][1] ≈  v
        @test df1[1][2] ≈  t
    end

    @testset "second order" begin
        pdp = Clapeyron.p∂p∂V(model,V,T,z)
        @test pdp[1] ≈ -v
        @test pdp[2] ≈ -vv
        h = Clapeyron.f_hess(model,V,T,z)
        @test all(h .≈ d2f)
        ddf = Clapeyron.∂2f(model,V,T,z)
        @test ddf[3] ≈ f
        @test all(ddf[2] .≈ df)
        @test all(ddf[1] .≈ d2f)
    end

    @testset "third order" begin
        dddf = Clapeyron.p∂p∂2p(model,V,T,z)
        @test dddf[2] ≈ -vv
        @test dddf[3] ≈ -vvv

        ddp = Clapeyron.∂2p(model,V,T,z)
        @test ddp[3] ≈ p
        @test all(ddp[2] .≈ dp)
        @test all(ddp[1] .≈ d2p)
    end
    
    @testset "AD Error Checks" begin
        import ForwardDiff
        import Clapeyron: MultipleTagError, NestedADError, __gradients_for_root_finders
        f_test(x,tups) = tups[1]*tups[2]*x # just a generic function to check error handling
        tag1 = ForwardDiff.Tag{:tag1,Float64}
        tag2 = ForwardDiff.Tag{:tag2,Float64}
        Tdual1 = ForwardDiff.Dual{tag1,Float64,1}
        Tdual2 = ForwardDiff.Dual{tag2,Float64,1}
        parts = ForwardDiff.Partials{1,Float64}((1.0,))
        theta1,theta2 = 0.5,2.0;
        x = 0.0
        x_dual = ForwardDiff.Dual{tag1,Float64,1}(x,parts)
        tups_primal = (theta1,theta2)
        tups_Dual2 = (Tdual1(theta1,parts),Tdual2(theta2,parts))
        tups_nestedDual = (ForwardDiff.Dual{tag2,Tdual1,1}(Tdual1(theta1,parts)),theta2)
        # Test multiple tag error
        @test_throws MultipleTagError __gradients_for_root_finders(x,tups_Dual2,tups_primal,f_test)
        # Test nested dual error
        @test_throws NestedADError __gradients_for_root_finders(x,tups_nestedDual,tups_primal,f_test)
        # Test dual as x error 
        @test_throws ErrorException __gradients_for_root_finders(x_dual,(theta1,theta2),tups_primal,f_test)
    end
end


