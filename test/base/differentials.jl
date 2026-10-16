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

@testset verbose = true "Differentials" begin
    
    @testset "Helmholtz bulk differentials" begin
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
    end


    @testset "Activity model differentials" begin
        #acetone-cloroform
        model_ge = VanLaar_GE(-0.8643,-0.5899) 
        model_gamma = VanLaar_lngamma(-0.8643,-0.5899)

        T = 300.0
        A12 = -0.8643 + 1e-5/T
        A21 = -0.5899 + 2e-5/T
        RT = Clapeyron.Rgas()*T
        zz = [0.4,1.6]
        n1, n2 = zz[1], zz[2]
        ax  = A12*n1 + A21*n2
        lnγ_exact = [A12*(A21*n2/ax)^2,A21*(A12*n1/ax)^2]
        gE_exact = RT*(A12*A21*n1*n2) / ax

        #∂lnγ/∂n
        pref = 2*A12^2*A21^2/ax^3
        ∂lnγ∂n_exact = pref * [ -n2^2   n1*n2 ;
                        n1*n2  -n1^2 ]

        #∂lnγ/∂T
        dA12dT = -1e-5/T^2
        dA21dT = -2e-5/T^2

        ∂lnγ∂T_exact = [ n2^2*A21*(-dA12dT*n1*A12*A21 + dA12dT*n2*A21^2 + 2*dA21dT*n1*A12^2)/ax^3,
                n1^2*A12*( 2*dA12dT*n2*A21^2 + dA21dT*n1*A12^2 - dA21dT*n2*A12*A21)/ax^3 ]


        cache1 = Clapeyron.∂lnϕ_cache(model_gamma,0.0,300.0,[0.4,1.6],Val(false))
        cache2 = Clapeyron.∂lnϕ_cache(model_gamma,0.0,300.0,[0.4,1.6],Val(true))
        cache3 = zeros(Float64,2)
        cache4 = Clapeyron.∂lnϕ_cache(model_ge,0.0,300.0,[0.4,1.6],Val(false)) #mainly to test errors
    
        #G_E
        g_E1 = Clapeyron.excess_gibbs_free_energy(model_ge,0.0,300.0,[0.4,1.6])
        g_E2 = Clapeyron.excess_gibbs_free_energy(model_gamma,0.0,300.0,[0.4,1.6])
        @test g_E1 ≈ g_E2 rtol = 1E-6
        @test g_E1 ≈ gE_exact rtol = 1E-6

        #lnγ
        lnγ1 = Clapeyron.lnγ(model_ge,0.0,300.0,[0.4,1.6])
        lnγ2 = Clapeyron.lnγ(model_gamma,0.0,300.0,[0.4,1.6])

        lnγ1 = Clapeyron.lnγ(model_ge,0.0,300.0,[0.4,1.6],cache1) |> deepcopy
        lnγ2 = Clapeyron.lnγ(model_gamma,0.0,300.0,[0.4,1.6],cache1) |> deepcopy

        lnγ1 = Clapeyron.lnγ(model_ge,0.0,300.0,[0.4,1.6],cache2) |> deepcopy
        lnγ2 = Clapeyron.lnγ(model_gamma,0.0,300.0,[0.4,1.6],cache2) |> deepcopy
        @test lnγ1 ≈ lnγ2 rtol = 1E-6
        @test lnγ1 ≈ lnγ_exact rtol = 1E-6

        lnγ1 = Clapeyron.lnγ(model_ge,0.0,300.0,[0.4,1.6],cache3) |> deepcopy
        lnγ2 = Clapeyron.lnγ(model_gamma,0.0,300.0,[0.4,1.6],cache3) |> deepcopy
        @test lnγ1 ≈ lnγ2 rtol = 1E-6
        @test lnγ1 ≈ lnγ_exact rtol = 1E-6

        #∂lnγ∂n
        g_E1,lnγ1,∂lnγ∂n1 = Clapeyron.∂lnγ∂n(model_ge,0.0,300.0,[0.4,1.6])
        g_E2,lnγ2,∂lnγ∂n2 = Clapeyron.∂lnγ∂n(model_gamma,0.0,300.0,[0.4,1.6])
        @test lnγ1 ≈ lnγ2 rtol = 1E-6
        @test lnγ1 ≈ lnγ_exact rtol = 1E-6
        @test g_E1 ≈ g_E2 rtol = 1E-6
        @test g_E1 ≈ gE_exact rtol = 1E-6
        @test ∂lnγ∂n1 ≈ ∂lnγ∂n2 rtol = 1E-6
        @test ∂lnγ∂n1 ≈ ∂lnγ∂n_exact rtol = 1E-6

        g_E1,_lnγ1,_∂lnγ∂n1 = Clapeyron.∂lnγ∂n(model_ge,0.0,300.0,[0.4,1.6],cache1)
        lnγ1,∂lnγ∂n1 = deepcopy(_lnγ1),deepcopy(_∂lnγ∂n1)
        g_E2,lnγ2,∂lnγ∂n2 = Clapeyron.∂lnγ∂n(model_gamma,0.0,300.0,[0.4,1.6],cache1)
        @test lnγ1 ≈ lnγ2 rtol = 1E-6
        @test lnγ1 ≈ lnγ_exact rtol = 1E-6
        @test g_E1 ≈ g_E2 rtol = 1E-6
        @test g_E1 ≈ gE_exact rtol = 1E-6
        @test ∂lnγ∂n1 ≈ ∂lnγ∂n2 rtol = 1E-6
        @test ∂lnγ∂n1 ≈ ∂lnγ∂n_exact rtol = 1E-6

        g_E1,_lnγ1,_∂lnγ∂n1 = Clapeyron.∂lnγ∂n(model_ge,0.0,300.0,[0.4,1.6],cache2)
        lnγ1,∂lnγ∂n1 = deepcopy(_lnγ1),deepcopy(_∂lnγ∂n1)
        g_E2,lnγ2,∂lnγ∂n2 = Clapeyron.∂lnγ∂n(model_gamma,0.0,300.0,[0.4,1.6],cache2)
        @test lnγ1 ≈ lnγ2 rtol = 1E-6
        @test lnγ1 ≈ lnγ_exact rtol = 1E-6
        @test g_E1 ≈ g_E2 rtol = 1E-6
        @test g_E1 ≈ gE_exact rtol = 1E-6
        @test ∂lnγ∂n1 ≈ ∂lnγ∂n2 rtol = 1E-6
        @test ∂lnγ∂n1 ≈ ∂lnγ∂n_exact rtol = 1E-6

        #∂lnγ∂n∂T
        g_E1,lnγ1,∂lnγ∂n1,∂lnγ∂T1 = Clapeyron.∂lnγ∂n∂T(model_ge,0.0,300.0,[0.4,1.6])
        g_E2,lnγ2,∂lnγ∂n2,∂lnγ∂T2 = Clapeyron.∂lnγ∂n∂T(model_gamma,0.0,300.0,[0.4,1.6])
        @test lnγ1 ≈ lnγ2 rtol = 1E-6
        @test lnγ1 ≈ lnγ_exact rtol = 1E-6
        @test g_E1 ≈ g_E2 rtol = 1E-6
        @test g_E1 ≈ gE_exact rtol = 1E-6
        @test ∂lnγ∂n1 ≈ ∂lnγ∂n2 rtol = 1E-6
        @test ∂lnγ∂n1 ≈ ∂lnγ∂n_exact rtol = 1E-6
        @test ∂lnγ∂T1 ≈ ∂lnγ∂T2 rtol = 1E-6
        @test ∂lnγ∂T1 ≈ ∂lnγ∂T_exact rtol = 1E-6

        g_E1,_lnγ1,_∂lnγ∂n1,_∂lnγ∂T1 = Clapeyron.∂lnγ∂n∂T(model_ge,0.0,300.0,[0.4,1.6],cache1)
        lnγ1,∂lnγ∂n1,∂lnγ∂T1 = deepcopy(_lnγ1),deepcopy(_∂lnγ∂n1),deepcopy(_∂lnγ∂T1)
        g_E2,lnγ2,∂lnγ∂n2,∂lnγ∂T2 = Clapeyron.∂lnγ∂n∂T(model_gamma,0.0,300.0,[0.4,1.6],cache1)
        @test lnγ1 ≈ lnγ2 rtol = 1E-6
        @test lnγ1 ≈ lnγ_exact rtol = 1E-6
        @test g_E1 ≈ g_E2 rtol = 1E-6
        @test g_E1 ≈ gE_exact rtol = 1E-6
        @test ∂lnγ∂n1 ≈ ∂lnγ∂n2 rtol = 1E-6
        @test ∂lnγ∂n1 ≈ ∂lnγ∂n_exact rtol = 1E-6
        @test ∂lnγ∂T1 ≈ ∂lnγ∂T2 rtol = 1E-6
        @test ∂lnγ∂T1 ≈ ∂lnγ∂T_exact rtol = 1E-6

        g_E1,_lnγ1,_∂lnγ∂n1,_∂lnγ∂T1 = Clapeyron.∂lnγ∂n∂T(model_ge,0.0,300.0,[0.4,1.6],cache2)
        lnγ1,∂lnγ∂n1,∂lnγ∂T1 = deepcopy(_lnγ1),deepcopy(_∂lnγ∂n1),deepcopy(_∂lnγ∂T1)
        g_E2,lnγ2,∂lnγ∂n2,∂lnγ∂T2 = Clapeyron.∂lnγ∂n∂T(model_gamma,0.0,300.0,[0.4,1.6],cache2)
        @test lnγ1 ≈ lnγ2 rtol = 1E-6
        @test lnγ1 ≈ lnγ_exact rtol = 1E-6
        @test g_E1 ≈ g_E2 rtol = 1E-6
        @test g_E1 ≈ gE_exact rtol = 1E-6
        @test ∂lnγ∂n1 ≈ ∂lnγ∂n2 rtol = 1E-6
        @test ∂lnγ∂n1 ≈ ∂lnγ∂n_exact rtol = 1E-6
        @test ∂lnγ∂T1 ≈ ∂lnγ∂T2 rtol = 1E-6
        @test ∂lnγ∂T1 ≈ ∂lnγ∂T_exact rtol = 1E-6
        @test_throws DimensionMismatch Clapeyron.∂lnγ∂n∂T(model_gamma,0.0,300.0,[0.4,1.6],cache4)

        #∂lnγ∂T
        ∂lnγ∂T1 = Clapeyron.∂lnγ∂T(model_ge,0.0,300.0,[0.4,1.6])
        ∂lnγ∂T2 = Clapeyron.∂lnγ∂T(model_gamma,0.0,300.0,[0.4,1.6])
        @test ∂lnγ∂T1 ≈ ∂lnγ∂T2 rtol = 1E-6
        @test ∂lnγ∂T1 ≈ ∂lnγ∂T_exact rtol = 1E-6

        ∂lnγ∂T1 = Clapeyron.∂lnγ∂T(model_ge,0.0,300.0,[0.4,1.6],cache1) |> deepcopy
        ∂lnγ∂T2 = Clapeyron.∂lnγ∂T(model_gamma,0.0,300.0,[0.4,1.6],cache1) |> deepcopy
        @test ∂lnγ∂T1 ≈ ∂lnγ∂T2 rtol = 1E-6
        @test ∂lnγ∂T1 ≈ ∂lnγ∂T_exact rtol = 1E-6

        ∂lnγ∂T1 = Clapeyron.∂lnγ∂T(model_ge,0.0,300.0,[0.4,1.6],cache2) |> deepcopy
        ∂lnγ∂T2 = Clapeyron.∂lnγ∂T(model_gamma,0.0,300.0,[0.4,1.6],cache2) |> deepcopy
        @test ∂lnγ∂T1 ≈ ∂lnγ∂T2 rtol = 1E-6
        @test ∂lnγ∂T1 ≈ ∂lnγ∂T_exact rtol = 1E-6
        @test_throws DimensionMismatch Clapeyron.∂lnγ∂T(model_gamma,0.0,300.0,[0.4,1.6],cache4)

        zz2 = 2 .* zz
        for cache in (nothing,cache1,cache2,cache3)
            for model in (model_ge, model_gamma)
                # G_E doubles
                g_E_1 = Clapeyron.excess_gibbs_free_energy(model, 0.0, 300.0, zz)
                g_E_2 = Clapeyron.excess_gibbs_free_energy(model, 0.0, 300.0, zz2)
                @test g_E_2 ≈ 2 * g_E_1 rtol = 1e-6

                # lnγ is unchanged
                _ln_gamma_1 = Clapeyron.lnγ(model, 0.0, 300.0, zz, cache)
                ln_gamma_1 = copy(_ln_gamma_1)
                ln_gamma_2 = Clapeyron.lnγ(model, 0.0, 300.0, zz2, cache)
                @test ln_gamma_2 ≈ ln_gamma_1 rtol = 1e-6

                if !(cache isa AbstractVector)
                    # ∂lnγ/∂T is unchanged
                    _dln_gamma_dT_1 = Clapeyron.∂lnγ∂T(model, 0.0, 300.0, zz, cache)
                    dln_gamma_dT_1 = copy(_dln_gamma_dT_1)
                    dln_gamma_dT_2 = Clapeyron.∂lnγ∂T(model, 0.0, 300.0, zz2, cache)
                    @test dln_gamma_dT_2 ≈ dln_gamma_dT_1 rtol = 1e-6

                    # ∂lnγ/∂n halves
                    g_E_1, _ln_gamma_1, _dln_gamma_dn_1 = Clapeyron.∂lnγ∂n(model, 0.0, 300.0, zz, cache)
                    dln_gamma_dn_1 = copy(_dln_gamma_dn_1)
                    ln_gamma_1 = copy(_ln_gamma_1)
                    g_E_2, _, dln_gamma_dn_2 = Clapeyron.∂lnγ∂n(model, 0.0, 300.0, zz2, cache)
                    @test g_E_2 ≈ 2 * g_E_1 rtol = 1e-6
                    @test ln_gamma_2 ≈ ln_gamma_1 rtol = 1e-6
                    @test dln_gamma_dn_2 ≈ dln_gamma_dn_1 ./ 2 rtol = 1e-6

                    # ∂lnγ/∂n∂T combined tuple consistency
                    g_E_1b, _ln_gamma_1b, _dln_gamma_dn_1b, _dln_gamma_dT_1b = Clapeyron.∂lnγ∂n∂T(model, 0.0, 300.0, zz, cache)
                    ln_gamma_1b, dln_gamma_dn_1b, dln_gamma_dT_1b = copy(_ln_gamma_1b), copy(_dln_gamma_dn_1b), copy(_dln_gamma_dT_1b)
                    g_E_2b, ln_gamma_2b, dln_gamma_dn_2b, dln_gamma_dT_2b = Clapeyron.∂lnγ∂n∂T(model, 0.0, 300.0, zz2, cache)

                    @test g_E_2b ≈ 2 * g_E_1b rtol = 1e-6
                    @test ln_gamma_2b ≈ ln_gamma_1b rtol = 1e-6
                    @test dln_gamma_dn_2b ≈ dln_gamma_dn_1b ./ 2 rtol = 1e-6
                    @test dln_gamma_dT_2b ≈ dln_gamma_dT_1b rtol = 1e-6
                end
            end
        end
    end
end
