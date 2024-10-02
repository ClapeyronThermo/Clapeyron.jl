using Clapeyron, Test

@printline

struct TestModel <: EoSModel end
Clapeyron.idealmodel(::TestModel) = BasicIdeal()
function Clapeyron.eos(model::TestModel,V,T,z)
    z_part = 1*z[1] + 2*z[2]+3*z[3]
    v_part = log(V)
    T_part = T^6
    vt_part = V*T*V*T
    return vt_part+T_part+v_part+z_part
end

function Clapeyron.∂f∂V(model::TestModel,V,T,z)
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
        dddf = Clapeyron.∂²³f(model,V,T,z)
        @test dddf[1] ≈ -vv
        @test dddf[2] ≈ -vvv

        ddp = Clapeyron.∂2p(model,V,T,z)
        @test ddp[3] ≈ p
        @test all(ddp[2] .≈ dp)
        @test all(ddp[1] .≈ d2p)
    end
end


