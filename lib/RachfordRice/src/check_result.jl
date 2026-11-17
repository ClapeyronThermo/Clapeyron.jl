module RachfordRiceConvergence

function log10_severity(residual, tolerance)
    iszero(residual) && return log10(zero(residual))
    return log10(residual / tolerance)
end

# ------------------------------------------------------------------
# Convergence checks
# ------------------------------------------------------------------
function composition_convergence_test(x,nc = length(x))
    residual         = abs(1 - sum(x))
    machine_epsilon  = eps(typeof(residual))
    absolute_tol     = 5 * machine_epsilon
    tolerance        = absolute_tol + nc * machine_epsilon
    severity         = log10_severity(residual, tolerance)
    return severity, severity <= 0
end

function phase_fraction_convergence_test(βV, βL)
    residual         = abs(βV + βL - 1) / (abs(βV) + abs(βL) + 1)
    machine_epsilon  = eps(typeof(residual))
    absolute_tol     = 5 * machine_epsilon
    severity         = log10_severity(residual, absolute_tol)
    return severity, severity <= 0
end

phase_fraction_convergence_test(βV) = phase_fraction_convergence_test(βV, 1 - βV)

function material_balance_convergence_test(βV, y, βL, x, z)
    w = z ./ sum(z)
    residual = maximum(
        abs.(βV .* y .+  βL .* x .- w) ./
        (abs.(βV .* y) .+ abs.(βL .* x) .+ w)
    )
    machine_epsilon  = eps(typeof(residual))
    absolute_tol     = 5 * machine_epsilon
    severity         = log10_severity(residual, absolute_tol)
    return severity, severity <= 0
end

function k_value_convergence_test(y, x, K)
    residual = maximum(
        abs.(y .- K .* x) ./
        (abs.(y) .+ abs.(K .* x))
    )
    machine_epsilon  = eps(typeof(residual))
    absolute_tol     = 5 * machine_epsilon
    severity         = log10_severity(residual, absolute_tol)
    return severity, severity <= 0
end

# ------------------------------------------------------------------
# Phase-boundary check
# ------------------------------------------------------------------
function within_phase_bounds(βV, K)
    Kmin,Kmax = extrema(K)
    min_βV = 1 / (1 - Kmax)
    max_βV = 1 / (1 - Kmin)
    return min_βV < βV < max_βV
end

# ------------------------------------------------------------------
# Main routine
# ------------------------------------------------------------------
function check_convergence(y, x, βV, βL, z, K)
    nc = length(K)
    y_sev, y_ok = composition_convergence_test(y, nc)
    x_sev, x_ok = composition_convergence_test(x, nc)

    fraction_sev, fraction_ok = phase_fraction_convergence_test(βV, βL)

    matbal_sev, matbal_ok =
        material_balance_convergence_test(βV, y, βL, x, z)

    k_sev, k_ok = k_value_convergence_test(y, x, K)

    bounds_ok = within_phase_bounds(βV, K)

    all_converged =
        y_ok && x_ok && fraction_ok && matbal_ok && k_ok && bounds_ok

    severities = (
        vapor_composition_severity   = y_sev,
        liquid_composition_severity  = x_sev,
        phase_fraction_severity      = fraction_sev,
        material_balance_severity    = matbal_sev,
        k_value_severity             = k_sev,
        all_converged                = all_converged,
    )

    convergences = (vapor_composition_converged  = y_ok,
        liquid_composition_converged = x_ok,
        phase_fraction_converged     = fraction_ok,
        material_balance_converged   = matbal_ok,
        k_value_converged            = k_ok,
        within_bounds                = bounds_ok
    )
    return all_converged, severities, convergences
end

end #module