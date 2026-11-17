const Conv = RachfordRice.RachfordRiceConvergence
using RachfordRice: residual, residual_first_derivative,
                    residual_first_second_derivative,
                    residual_first_second_derivative_error,
                    rr_limits, rr_status

function test_cases(CASES)
    failures = []
    for case in CASES
        K,z = case.K,case.z
        RR = RachfordRiceResidual(K,z)
        result = RachfordRice.rachfordrice(RR,K_tol = 0.0,negative_flash = true)
        L,V = result.L, result.V
        x = liquid_compositions(RR,result)
        y = vapour_compositions(RR,result)
        all_converged,severities, convergences = Conv.check_convergence(y, x, V, L, z, K)
        if !all_converged
            push!(failures,(RR, result, severities, convergences))
        end
    end
    if length(failures) == 0
        @test length(failures) == 0
    else
        for failure in failures
            (_, result, severities, convergences) = failure
            @test result.status != RRFailure
            @test severities.vapor_composition_severity < 5
            @test severities.liquid_composition_severity < 5
            @test severities.phase_fraction_severity <= 0
            @test severities.material_balance_converged <= 0
            @test severities.k_value_severity <= 0
            @test convergences.within_bounds
        end
    end
    return failures
end