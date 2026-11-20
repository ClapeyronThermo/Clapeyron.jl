using Downloads
using CSV
using RachfordRice
const Conv = RachfordRice.RachfordRiceConvergence
using RachfordRice: residual, residual_first_derivative,
                    residual_first_second_derivative,
                    residual_first_second_derivative_error,
                    rr_limits, rr_status