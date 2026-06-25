
function bubble_point(model::EoSModel,T,x)
    return 1,1
end

function ff_582 end
bulk_p(model, T, v) = pressure(model, v, T)
two_props(model, T, v) = (pressure(model, v, T), Clapeyron.molecular_weight(model) / 1e-3)
mse_loss(y_calc, y_exp) = abs2(y_calc - y_exp)

# Tiny in-memory dataset (single output)
function make_pv_data(model, Ts, vs)
    ps = [pressure(model, v, T) for (T, v) in zip(Ts, vs)]
    (T = Ts, v = vs, out_p = ps)
end


@testset "Estimation framework" begin

    @testset "__mse default loss" begin
        @test Clapeyron.__mse(1.1, 1.0) ≈ abs2((1.1 - 1.0) / 1.0)
        @test Clapeyron.__mse(1.0, 1.0) == 0.0
        @test Clapeyron.__mse(0.5, 1.0) ≈ abs2((0.5 - 1.0) / 1.0)
    end

    # ------------------------------------------------------------------
    @testset "EstimationModel" begin
        model = PCSAFT("methane")

        toestimate = [
            Dict(:param => :epsilon, :lower => 100., :upper => 300., :guess => 250.),
            Dict(:param => :sigma,   :factor => 1e-10, :lower => 3.2, :upper => 4.0, :guess => 3.7),
            Dict(:param => :segment, :lower => 0.9, :upper => 1.1, :guess => 1.),
        ]

        est = EstimationModel(model, toestimate)
        test_repr(est)
        test_repr(est.toestimate)
        #construction
        @test est isa EstimationUtils.AbstractEstimationModel
        @test EstimationUtils.get_model(est) === est.model

        #parameter_length
        @test EstimationUtils.parameter_length(est) == 3

        #bounds and initial_guess
        lb = EstimationUtils.lower_bounds(est)
        ub = EstimationUtils.upper_bounds(est)
        x0 = EstimationUtils.initial_guess(est)

        @test lb == [100., 3.2, 0.9]
        @test ub == [300., 4.0, 1.1]
        @test x0 == [250., 3.7, 1.0]
        @test length(lb) == length(ub) == length(x0) == 3

        #get_eos_parameters / set_eos_parameters! round-trip
        Θ = EstimationUtils.get_eos_parameters(est)
        @test length(Θ) == 3

        # Modify the epsilon value (index 1, no factor)
        Θ_new = copy(Θ)
        Θ_new[1] = 200.0
        EstimationUtils.set_eos_parameters!(est, Θ_new)
        Θ_back = EstimationUtils.get_eos_parameters(est)
        @test Θ_back[1] ≈ 200.0


        # Modify the sigma value (index 2, factor = 1e-10)
        # The value exposed to the optimizer is in Å (÷ 1e-10)
        Θ_new2 = copy(Θ_back)
        Θ_new2[2] = 3.5   # angstroms in optimizer space
        EstimationUtils.set_eos_parameters!(est, Θ_new2)
        Θ_back2 = EstimationUtils.get_eos_parameters(est)
        @test Θ_back2[2] ≈ 3.5   # round-trip through factor
        # The underlying parameter is stored in SI (meters)
        @test est.model.params.sigma[1] ≈ 3.5e-10

        #factor: SI storage vs optimizer space
        # sigma is stored in meters; the optimizer sees it in units of 1e-10 m
        Θ = EstimationUtils.get_eos_parameters(est)
        # Whatever the current sigma is in meters, optimizer sees it divided by 1e-10
        sigma_si = est.model.params.sigma[1]
        @test Θ[2] ≈ sigma_si / 1e-10 rtol = 1e-8

        #indexing interface
        est[1] = 150.0
        @test est[1] ≈ 150.0
        est[2] = 3.8
        @test est[2] ≈ 3.8
        est[3] = 1.05
        @test est[3] ≈ 1.05

        #broadcasting (copyto!)
        Θ_bc = [160.0, 3.6, 0.95]
        est .= Θ_bc
        Θ_read = EstimationUtils.get_eos_parameters(est)
        @test Θ_read ≈ Θ_bc

        #symbol indexing
        est[:epsilon] = 312
        m = EstimationUtils.get_model(est)
        @test m.params.epsilon[1] ≈ 312
        
        est[[:sigma,:epsilon]] = (2.5,313)
        @test est[:sigma] ≈ [2.5]
        @test est[:epsilon] ≈ [313]
        
        #set_model / get_model defaults
        
        model2 = PCSAFT("methane")
        est2 = EstimationUtils.set_model(est, model2)
        @test EstimationUtils.get_model(est2) === model2
        @test EstimationUtils.get_model(est) !== model2

        #default bounds when not specified
        toestimate_nobound = [Dict(:param => :epsilon)]
        est_nb = EstimationModel(model, toestimate_nobound)
        lb = EstimationUtils.lower_bounds(est_nb)
        ub = EstimationUtils.upper_bounds(est_nb)
        @test all(lb .== -Inf)
        @test all(ub .== +Inf)

        #initial_guess falls back to model value when :guess omitted
        toestimate_noguess = [Dict(:param => :epsilon, :lower => 100., :upper => 300.)]
        est_ng = EstimationModel(model, toestimate_noguess)
        x0 = EstimationUtils.initial_guess(est_ng)
        # Should equal the epsilon currently stored in the model
        @test x0[1] ≈ est_ng.model.params.epsilon[1]
    end

    @testset "EstimationModel – symmetric" begin
        model = NRTL(["ethanol", "methanol"],puremodel = BasicIdeal)

        toestimate = [
            Dict(
                :param   => :a,
                :indices => (1, 2),
                :symmetric => true,
                :lower   => 100,
                :upper   =>  500,
                :guess   =>  200,
            )
        ]

        toestimate2 = [
            Dict(
                :param   => :a,
                :indices => (1, 2),
                :symmetric => false,
                :lower   => 100,
                :upper   =>  500,
                :guess   =>  200,
            )
        ]
        sym_est = EstimationModel(model, toestimate)
        asym_est = EstimationModel(model, toestimate2)
        #set then get round-trip for off-diagonal pair
        EstimationUtils.set_eos_parameters!(sym_est, [250.0])
        @test model.params.a.values[1,2] == 250.0
        @test model.params.a.values[2,1] == 250.0

        EstimationUtils.set_eos_parameters!(asym_est, [350.0])
        @test model.params.a.values[1,2] == 350.0
        @test model.params.a.values[2,1] == 250.0
    end

    @testset "EstimationModel – cross_assoc" begin
        model = SAFTVRMie(["ethanol", "water"],assoc_options = :cr1)
        ap = model.params.epsilon_assoc   # AssocParam
        
        # Snapshot values we must not touch
        val1_before = ap.values.values[1]

        toestimate = [
            Dict(
                :param      => :epsilon_assoc,
                :indices    => 2,
                :cross_assoc => true,
                :lower      => 2000.,
                :upper      => 4000.,
                :guess      => 3500.,
            )
        ]

        toestimate2 = [
            Dict(
                :param      => :epsilon_assoc,
                :indices    => 2,
                :cross_assoc => false,
                :lower      => 2000.,
                :upper      => 4000.,
                :guess      => 3500.,
            )
        ]

        assoc_est = EstimationModel(model, toestimate)
        assoc_est2 = EstimationModel(model, toestimate2)
        #parameter_length is 1
        @test EstimationUtils.parameter_length(assoc_est) == 1

        #set_eos_parameters! writes the primary value
        new_val = 3200.0
        EstimationUtils.set_eos_parameters!(assoc_est, [new_val])
        @test ap.values.values[2] ≈ new_val

        #cross_assoc mirrors value to the symmetric inner entry
        # Retrieve the outer/inner indices that define which pair matrix is touched.
        k     = 2
        ij    = ap.values.outer_indices[k]
        ab    = ap.values.inner_indices[k]
        i, j  = ij
        a, b  = ab

        ij_mat = ap.values[i, j]
        new_val = 3100.0
        EstimationUtils.set_eos_parameters!(assoc_est, [new_val])
        #with symmetric cross-assoc, Both [a,b] and [b,a] must equal new_val
        @test ij_mat[a, b] ≈ new_val
        @test ij_mat[b, a] ≈ new_val

        EstimationUtils.set_eos_parameters!(assoc_est2, [2*new_val])
        #without symmetric cross-assoc, ij[a,b] == ij[b,a]/2
        @test ij_mat[a, b] ≈ 2*new_val
        @test ij_mat[b, a] ≈ new_val

        #other association values are not disturbed
        @test ap.values.values[1] ≈ val1_before
    end

    @testset "EstimationModel - multiple-value indices" begin
        model_multiple = PCSAFT(["methane","ethane"])

        toestimate_multiple = [
            Dict(
                :param   => :epsilon,
                :indices => :pures,       # one Θ entry per component: ε₁, ε₂
            ),
            Dict(
                :param     => :epsilon,
                :indices   => :unlike,    # one Θ entry for the cross term ε₁₂
                :symmetric => true,       # ε₁₂ == ε₂₁ (default, shown for clarity)
            ),
            Dict(
                :param   => :sigma,
                :factor  => 1e-10,
                :indices => :pures,       # σ₁, σ₂
            ),
        ]
        # Θ has 5 elements: [ε₁, ε₂, ε₁₂, σ₁, σ₂]
        est_model = EstimationModel(model_multiple, toestimate_multiple)
        @test EstimationUtils.parameter_length(est_model) == 5
        est_model[3] = 170.0
        @test model_multiple.params.epsilon[1,2] == 170.0
    end

    @testset "EstimationData – from NamedTuple" begin
        model = cPR("methane")
        Ts = [300., 350., 400.]
        vs = [1e-3, 1.2e-3, 1.4e-3]
        table = make_pv_data(model, Ts, vs)
        est_data = EstimationData(table, bulk_p, mse_loss)
        test_repr(est_data)
        #construction
        @test est_data isa EstimationUtils.AbstractEstimationLoss
        @test length(est_data.inputs) == 3
        @test length(est_data.outputs) == 3
        @test est_data.inputs_name == [:T, :v]
        @test est_data.outputs_name == [:out_p]
        @test est_data.valid == [true,true,true]

        #zero-inputs
        est_data_0 = EstimationData((Tc = [1.0],Vc = [1.0]),bulk_p,mse_loss)
        @test est_data_0.valid == [true]

        #missing values
        est_data_missing = EstimationData((a = [1.0,missing,1.0],b = [1.0,1.0,missing],out_c = [1.0,1.0,1.0]),bulk_p,mse_loss)
        @test est_data_missing.valid == [true,false,false]
        @test isequal(est_data_missing.inputs,[(1.0,1.0),(NaN,1.0),(1.0,NaN)])
        @test est_data_missing.inputs_ismissingvalues == [(false,false),(true,false),(false,true)]

        #objective_function – zero loss at reference model
        # model produced the data so the loss must be zero
        F = EstimationUtils.objective_function(est_data, model)
        @test F ≈ 0.0 atol = 1e-10

        #objective_function – non-zero loss on different model
        model2 = PCSAFT("methane")
        F_model2 = EstimationUtils.objective_function(est_data, model2)
        @test F_model2 > 0.0
        @test isfinite(F_model2)

        #normalize flag
        # With normalize = true (default) the result is divided by npoints
        est_norm   = EstimationData(table, bulk_p, mse_loss; normalize = true)
        est_nonorm = EstimationData(table, bulk_p, mse_loss; normalize = false)
        model2 = PCSAFT("methane")
        F_norm   = EstimationUtils.objective_function(est_norm,   model2)
        F_nonorm = EstimationUtils.objective_function(est_nonorm, model2)
        @test F_nonorm ≈ F_norm * length(est_norm.inputs) rtol = 1e-8

        #data_weights – scalar multiplication
        w = 2.0
        data_w = fill(w, length(Ts))
        est_w  = EstimationData(table, bulk_p, mse_loss; data_weights = data_w, normalize = false)
        est_uw = EstimationData(table, bulk_p, mse_loss; normalize = false)
        model2 = PCSAFT("methane")
        Fw  = EstimationUtils.objective_function(est_w,  model2)
        Fuw = EstimationUtils.objective_function(est_uw, model2)
        @test Fw ≈ w * Fuw rtol = 1e-8
    end

    @testset "EstimationData – multiple outputs" begin
        model = cPR("methane")
        Ts = [300., 350., 400.]
        vs = [1e-3, 1.2e-3, 1.4e-3]

        ps = [pressure(model, v, T) for (T, v) in zip(Ts, vs)]
        mws = fill(Clapeyron.molecular_weight(model) / 1e-3, 3)   # constant, known value

        table = (T = Ts, v = vs, out_p = ps, out_mw = mws)
        est_data = EstimationData(table, two_props, mse_loss)

        @test length(est_data.outputs_name) == 2
        @test :out_p  in est_data.outputs_name
        @test :out_mw in est_data.outputs_name

        F = EstimationUtils.objective_function(est_data, model)
        @test F ≈ 0.0 atol = 1e-10

        est_w1 = EstimationData(table, two_props, mse_loss; output_weights = (1.0, 1.0), normalize = false)
        est_w2 = EstimationData(table, two_props, mse_loss; output_weights = (2.0, 2.0), normalize = false)
        model2 = PCSAFT("methane")
        F1 = EstimationUtils.objective_function(est_w1, model2)
        F2 = EstimationUtils.objective_function(est_w2, model2)
        @test F2 ≈ 2.0 * F1 rtol = 1e-8
    end

    @testset "EstimationData – CSV option parsing" begin
        # Check that the JSON header format is parsed correctly (re: issue #582)
        csv_json = """Clapeyron Estimator,,,
        {"method" : "bulk_p","species" : ["carbon dioxide","ethane"]},,
        T,v,out_p"""

        ed = Clapeyron.EstimationData(csv_json)
        @test ed.species == ["carbon dioxide", "ethane"]

        # Check vec-style header format
        csv_vec = """Clapeyron Estimator,,,
        [method = bulk_p,species = "carbon dioxide" "ethane",]
        T,v,out_p"""

        ed2 = Clapeyron.EstimationData([csv_vec])
        @test ed2.species == ["carbon dioxide", "ethane"]
    end

    @testset "EstimationData – error column recognition" begin
        # Columns ending in _error_abs, _error_rel, _error_std should NOT be
        # treated as inputs or outputs.
        model = cPR("methane")
        Ts = [300., 350.]
        vs = [1e-3, 1.1e-3]
        ps = [pressure(model, v, T) for (T, v) in zip(Ts, vs)]
        p_err = [1e2, 1e2]

        table = (T = Ts, v = vs, out_p = ps, out_p_error_abs = p_err)
        est_data = EstimationData(table, bulk_p, mse_loss)
        @test :out_p_error_abs ∉ est_data.outputs_name
        @test :out_p           in  est_data.outputs_name
        # The error for the single output should be filled in
        @test all(e[1] ≈ 1e2 for e in est_data.outputs_error)
    end

    @testset "EstimationData - edge cases" begin
        csv_582 = """Clapeyron Estimator,,,
[method = sum,species=acetonitrile heptane palmitic acid],,
T,z1,z2,z3,out_x1,out_x2,out_x3,out_y1,out_y2,out_y3"""
        data = EstimationData(csv_582)
        model = MonomerIdeal(["acetonitrile","heptane","palmitic acid"])
        est_model = EstimationModel(model,[Dict(:param => :Mw,:indices => :all)])

        @test_throws "Try wrapping each species name in the list" EstimationProblem(est_model,[data])
    end

    @testset "EstimationProblem" begin
        model = PCSAFT("methane")
        toestimate = [
            Dict(:param => :epsilon, :lower => 100., :upper => 300., :guess => 250.),
            Dict(:param => :sigma,   :factor => 1e-10, :lower => 3.2, :upper => 4.0, :guess => 3.7),
            Dict(:param => :segment, :lower => 0.9, :upper => 1.1, :guess => 1.),
        ]
        est_model = EstimationModel(model, toestimate)

        ref_model = cPR("methane")
        Ts = [300., 350., 400.]
        vs = [1e-3, 1.2e-3, 1.4e-3]
        table = make_pv_data(ref_model, Ts, vs)
        est_data = EstimationData(table, bulk_p, mse_loss)

        prob = EstimationProblem(est_model, [est_data])

        #bounds and initial_guess are forwarded" begin
        @test EstimationUtils.lower_bounds(prob)  == EstimationUtils.lower_bounds(est_model)
        @test EstimationUtils.upper_bounds(prob)  == EstimationUtils.upper_bounds(est_model)
        @test EstimationUtils.initial_guess(prob) == EstimationUtils.initial_guess(est_model)

        #objective_function(prob, Θ)
        Θ0 = EstimationUtils.initial_guess(prob)
        F = EstimationUtils.objective_function(prob, Θ0)
        @test isfinite(F)
        @test F >= 0.0

        #objective_function(prob) returns a callable
        f = EstimationUtils.objective_function(prob)
        @test f isa Function
        Θ0 = EstimationUtils.initial_guess(prob)
        @test f(Θ0) ≈ EstimationUtils.objective_function(prob, Θ0)

        #multiple data objects: total loss is sum
        # Build a second dataset with the same data; total must be 2× single
        est_data2 = EstimationData(table, bulk_p, mse_loss)
        prob1 = EstimationProblem(est_model, [est_data])
        prob2 = EstimationProblem(est_model, [est_data, est_data2])
        Θ0 = EstimationUtils.initial_guess(prob1)
        F1 = EstimationUtils.objective_function(prob1, Θ0)
        F2 = EstimationUtils.objective_function(prob2, Θ0)
        @test F2 ≈ 2.0 * F1 rtol = 1e-8

        #parameter update is reflected in underlying model
        Θ_new = [200., 3.5, 1.0]
        EstimationUtils.set_eos_parameters!(prob.toestimate, Θ_new)
        Θ_back = EstimationUtils.get_eos_parameters(prob.toestimate)
        @test Θ_back ≈ Θ_new rtol = 1e-8
    end

    @testset "Custom AbstractEstimationModel – KijEstimationModel" begin
        # Minimal stub matching the interface shown in custom_estimation.md.
        # We use a binary PR model; cubic models have get_k / set_k!.

        struct KijEstimationModel{M<:Clapeyron.CubicModel} <: EstimationUtils.AbstractEstimationModel{M}
            model::M
        end

        function EstimationUtils.get_eos_parameters(est::KijEstimationModel)
            [Clapeyron.get_k(est.model)[1, 2]]
        end

        function EstimationUtils.set_eos_parameters!(est::KijEstimationModel, k12_vec)
            k12 = only(k12_vec)
            # build symmetric 2×2 matrix
            Clapeyron.set_k!(est.model, [0. k12; k12 0.])
        end

        EstimationUtils.initial_guess(est::KijEstimationModel) = [0.0]
        EstimationUtils.lower_bounds(est::KijEstimationModel)  = [-1.0]
        EstimationUtils.upper_bounds(est::KijEstimationModel)  = [ 1.0]

        model = PR(["carbon dioxide", "methane"])
        est   = KijEstimationModel(model)

        #get_model default
        @test EstimationUtils.get_model(est) === model

        #initial_guess
        @test EstimationUtils.initial_guess(est) == [0.0]

        #set_eos_parameters! / get_eos_parameters round-trip
        EstimationUtils.set_eos_parameters!(est, [0.1])
        @test Clapeyron.get_k(model)[1, 2] ≈ 0.1
        Θ = EstimationUtils.get_eos_parameters(est)
        @test Θ ≈ [0.1]

        #bounds
        @test EstimationUtils.lower_bounds(est) == [-1.0]
        @test EstimationUtils.upper_bounds(est) == [ 1.0]

        #parameter_length falls back to length of get_eos_parameters
        EstimationUtils.set_eos_parameters!(est, [0.0])  # reset
        @test EstimationUtils.parameter_length(est) == 1

        #set_model / get_model
        model2 = PR(["carbon dioxide", "methane"])
        est2   = EstimationUtils.set_model(est, model2)
        @test EstimationUtils.get_model(est2) === model2
        @test EstimationUtils.get_model(est)  === model
    end

    @testset "Custom AbstractEstimationLoss – bulk property version" begin
        # Stub matching the AbstractEstimationLoss interface from custom_estimation.md.
        # Instead of saturation pressure (phase equilibrium), we use pressure() directly.

        struct BulkPressureLoss <: EstimationUtils.AbstractEstimationLoss
            T_exp::Vector{Float64}
            v_exp::Vector{Float64}
            p_ref::Vector{Float64}   # reference pressures pre-computed
        end

        function BulkPressureLoss(ref_model::Clapeyron.EoSModel, Ts, vs)
            p_ref = [pressure(ref_model, v, T) for (T, v) in zip(Ts, vs)]
            BulkPressureLoss(Ts, vs, p_ref)
        end

        function EstimationUtils.objective_function(loss::BulkPressureLoss, model)
            F     = zero(eltype(model))
            valid = 0
            for (T, v, p_ref) in zip(loss.T_exp, loss.v_exp, loss.p_ref)
                p_calc = pressure(model, v, T)
                if isfinite(p_calc) && isfinite(p_ref)
                    F     += abs2((p_calc - p_ref) / p_ref)
                    valid += 1
                end
            end
            return valid > 0 ? F / valid : Inf * one(typeof(F))
        end

        ref_model = cPR("methane")
        Ts = [300., 350., 400., 450.]
        vs = [1e-3, 1.1e-3, 1.2e-3, 1.3e-3]
        loss = BulkPressureLoss(ref_model, Ts, vs)
        test_repr(loss)
        #zero loss at reference model
        F = EstimationUtils.objective_function(loss, ref_model)
        @test F ≈ 0.0 atol = 1e-10

        #positive loss on different model
        candidate = PCSAFT("methane")
        F = EstimationUtils.objective_function(loss, candidate)
        @test F > 0.0
        @test isfinite(F)

        #composes into EstimationProblem
        model = cPR("methane")
        toestimate = [Dict(:param => :a, :lower => 0.1, :upper => 1.0)]
        est_model = EstimationModel(model, toestimate)
        prob = EstimationProblem(est_model, [loss])
        Θ0 = EstimationUtils.initial_guess(prob)
        F = EstimationUtils.objective_function(prob, Θ0)
        @test isfinite(F)
        @test F >= 0.0
    end

    @testset "Composed problem – multiple heterogeneous loss terms" begin
        # Mirrors the pattern shown at the end of custom_estimation.md:
        # one standard EstimationData + one custom AbstractEstimationLoss,
        # both minimised simultaneously.

        model = PCSAFT("methane")
        toestimate = [
            Dict(:param => :epsilon, :lower => 100., :upper => 400., :guess => 200.),
        ]
        est_model = EstimationModel(model, toestimate)

        ref = cPR("methane")
        Ts = [300., 350.]
        vs = [1e-3, 1.1e-3]
        table = make_pv_data(ref, Ts, vs)

        loss_A = EstimationData(table, bulk_p, mse_loss)

        struct RefPressureLoss2 <: EstimationUtils.AbstractEstimationLoss
            T::Vector{Float64}
            v::Vector{Float64}
            pref::Vector{Float64}
        end

        function EstimationUtils.objective_function(l::RefPressureLoss2, m)
            F = zero(eltype(m))
            for (T, v, p0) in zip(l.T, l.v, l.pref)
                p = pressure(m, v, T)
                F += abs2((p - p0)/p0)
            end
            F / length(l.T)
        end

        pref = [pressure(ref, v, T) for (T, v) in zip(Ts, vs)]
        loss_B = RefPressureLoss2(Ts, vs, pref)

        prob_single = EstimationProblem(est_model, [loss_A])
        prob_double = EstimationProblem(est_model, [loss_A, loss_B])
        test_repr(prob_single)
        test_repr(prob_double)
        Θ0 = EstimationUtils.initial_guess(prob_single)

        FA = EstimationUtils.objective_function(prob_single, Θ0)
        FAB = EstimationUtils.objective_function(prob_double, Θ0)

        # loss_B at the same Θ0
        EstimationUtils.set_eos_parameters!(prob_single.toestimate, Θ0)
        FB = EstimationUtils.objective_function(loss_B, prob_single.model)

        @test FAB ≈ FA + FB rtol = 1e-8
    end

    @testset "Estimation - SAFTgammaMie" begin
        model = SAFTgammaMie(["ethanol","water"],epsilon_mixing = :hudsen_mccoubrey)

        toestimate = [Dict(
            :param => :epsilon,
            :indices => (1,1),
            :recombine => true,
            :lower => 100.,
            :upper => 300.,
            :guess => 250.
        ),
        Dict(
            :param => :epsilon,
            :indices => (2,3),
            :symmetric => false,
            :lower => 100.,
            :upper => 300.,
            :guess => 250.
        ),
        Dict(
            :param => :sigma,
            :indices => (1,1),
            :factor => 1e-10,
            :lower => 3.2,
            :upper => 4.0,
            :guess => 3.7
            ),
        Dict(
            :param => :epsilon_assoc,
            :indices => 2,
            :cross_assoc => true,
            :lower => 2000.,
            :upper => 4000.,
            :guess => 3500.
            )]

        estimator,objective,initial,upper,lower = Estimation(model,toestimate,["../examples/data/bubble_point.csv"],[:vrmodel])

        model2 = return_model(estimator,model,initial)
        @test model2.params.epsilon[1,1] == initial[1] # Test that the parameter was correctly updated
        @test model2.params.epsilon[1,2] ≈ 251.57862124740765 rtol = 1e-6 # Test that the combining rule was used
        @test model2.params.epsilon[1,2] == model2.params.epsilon[2,1] # Test that the unlike parameters remain symmetric

        @test model2.params.epsilon[2,3] == initial[2] # Test that the parameter was updated
        @test model2.params.epsilon[2,3] != model2.params.epsilon[3,2] # Test that the parameter is no longer symmetric

        @test model2.params.sigma[1,1] == initial[3]*1e-10 # Test that the factor was used to update the parameters

        @test model2.params.epsilon[2,3] == initial[2] # Test that the parameter was updated
        @test model2.params.epsilon[2,3] != model2.params.epsilon[3,2] # Test that the parameter is no longer symmetric

        @test model2.params.epsilon_assoc.values.values[2] == initial[4] # Test that the association parameter was updated
        @test model2.params.epsilon_assoc.values.values[2] == model2.params.epsilon_assoc.values.values[3] # Test that the cross-association parameter was updated

        @test objective(initial) ≈ 2.4184631612655836 rtol = 1e-6

        estimator,objective,initial,upper,lower = Estimation(model,toestimate,[(2.,"../examples/data/bubble_point.csv"),(1.,"../examples/data/bubble_point.csv")],[:vrmodel])

        @test objective(initial) ≈ 7.255389483796751 rtol = 1e-6
    end
end #@testset "Estimation Framework"


@testset "promote_model" begin
    @testset "#365" begin
        #error found during #365
        model = SAFTgammaMie(["ethanol","water"],epsilon_mixing = :hudsen_mccoubrey)
        modelvec = Clapeyron.EoSVectorParam(model)
        @test Clapeyron.promote_model(BigFloat,modelvec) isa Clapeyron.EoSVectorParam
    end

    @testset "#366" begin
        #=
        #366
        incorrect conversion of MixedGCSegmentParam.
        =#
        model = SAFTgammaMie(["ethanol","water"],epsilon_mixing = :hudsen_mccoubrey)
        bigfloat_model = Clapeyron.promote_model(BigFloat,model)
        @test model.params.mixed_segment.values.v ≈ bigfloat_model.params.mixed_segment.values.v

        model2 = SAFTgammaMie(["octane"])
        k1 = copy(model2.params.mixed_segment.values.v)
        Clapeyron.recombine!(model2)
        k2 = model2.params.mixed_segment.values.v
        @test k1[1] == k2[1] # Test that the first value is unchanged
    end
end

#Estimation(model,deepcopy(model),Clapeyron.ToEstimate(toestimate_582),[:puremodel],Clapeyron.__mse)