const CLAPEYRON_LIB_PATH = normpath(Base.pkgdir(Clapeyron), "lib")

const RR_TEST_PATH = joinpath(CLAPEYRON_LIB_PATH, "RachfordRice", "test")

using Clapeyron.RachfordRice

const Conv = Clapeyron.RachfordRice.RachfordRiceConvergence
using Clapeyron.RachfordRice: residual, residual_first_derivative, residual_first_second_derivative, residual_first_second_derivative_error, rr_limits, rr_status

include(joinpath(RR_TEST_PATH, "helpers.jl"))
include(joinpath(RR_TEST_PATH, "whitson_rr_contest.jl"))
include(joinpath(RR_TEST_PATH, "extra_cases.jl"))
include(joinpath(RR_TEST_PATH, "api.jl"))
