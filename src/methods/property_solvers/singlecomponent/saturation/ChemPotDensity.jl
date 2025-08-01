function ∂Helmholtz(model::EoSModel, ρ, T, z=[1.0])
    # Auxiliar functions to compute Helmholtz energy and its first density derivative
    f(dρ) = eos(model, 1. / dρ, T, z)
    A, ∂A, = Solvers.f∂f(f, ρ)
    return A, ∂A
end

function ∂2Helmholtz(model::EoSModel, ρ, T, z=[1.0])
    # Auxiliar functions to compute Helmholtz energy and its first and second density derivatives
    f(dρ) = eos(model, 1. / dρ, T, z)
    A, ∂A, ∂2A = Solvers.f∂f∂2f(f, ρ)
    return A, ∂A, ∂2A
end

function fobj_psat!(model::EoSModel, T)
    
    # Objetive function to solve saturation pressure using the densities as iterable variables
    # T = Saturation Temperature
    # ρ = initial guess for the phases densities = [ρ liquid, ρ vapor]
    # F = vector for objective function
    # J = matrix for objective function jacobian

    ps,μs = equilibria_scale(model)
    function f!(F,x)
        ρ_liq, ρ_vap = x
        A_liq, ∂A_liq = ∂Helmholtz(model, ρ_liq, T)
        A_vap, ∂A_vap = ∂Helmholtz(model, ρ_vap, T)
        P_liq = ∂A_liq*ρ_liq^2
        μ_liq = A_liq + ρ_liq * ∂A_liq
        P_vap = ∂A_vap*ρ_vap^2
        μ_vap = A_vap + ρ_vap * ∂A_vap
        F[1] = (μ_liq - μ_vap)*μs
        F[2] = (P_liq - P_vap)*ps
        return F
    end

    function fj!(F,J,x)
        ρ_liq, ρ_vap = x
        A_liq, ∂A_liq, ∂2A_liq = ∂2Helmholtz(model, ρ_liq, T)
        ∂P_liq = 2. * ρ_liq * ∂A_liq + ρ_liq^2 * ∂2A_liq
        ∂μ_liq = ρ_liq*∂2A_liq + 2*∂A_liq
        A_vap, ∂A_vap, ∂2A_vap = ∂2Helmholtz(model, ρ_vap, T)
        ∂P_vap = 2. * ρ_vap * ∂A_vap + ρ_vap^2 * ∂2A_vap
        ∂μ_vap = ρ_vap*∂2A_vap + 2*∂A_vap
        J[1, 1] = ∂μ_liq*μs
        J[1, 2] = -∂μ_vap*μs
        J[2, 1] = ∂P_liq*ps
        J[2, 2] = -∂P_vap*ps
        P_liq = ∂A_liq*ρ_liq^2
        μ_liq = A_liq + ρ_liq * ∂A_liq
        P_vap = ∂A_vap*ρ_vap^2
        μ_vap = A_vap + ρ_vap * ∂A_vap
        F[1] = (μ_liq - μ_vap)μs
        F[2] = (P_liq - P_vap)*ps
        return F,J
    end

    function j!(J,x)
        ρ_liq, ρ_vap = x
        A_liq, ∂A_liq, ∂2A_liq = ∂2Helmholtz(model, ρ_liq, T)
        ∂P_liq = 2. * ρ_liq * ∂A_liq + ρ_liq^2 * ∂2A_liq
        ∂μ_liq = ρ_liq*∂2A_liq + 2*∂A_liq
        A_vap, ∂A_vap, ∂2A_vap = ∂2Helmholtz(model, ρ_vap, T)
        ∂P_vap = 2. * ρ_vap * ∂A_vap + ρ_vap^2 * ∂2A_vap
        ∂μ_vap = ρ_vap*∂2A_vap + 2*∂A_vap
        J[1, 1] = ∂μ_liq*μs
        J[1, 2] = -∂μ_vap*μs
        J[2, 1] = ∂P_liq*ps
        J[2, 2] = -∂P_vap*ps
        return J
    end

    function jv!(x)
       return nothing
    end
    return NLSolvers.VectorObjective(f!,j!,fj!,jv!) |> NLSolvers.NEqProblem
end


struct ChemPotDensitySaturation{T,C} <: SaturationMethod
    vl::Union{Nothing,T}
    vv::Union{Nothing,T}
    crit::C
    f_limit::Float64
    atol::Float64
    rtol::Float64
    max_iters::Int
end

"""
    ChemPotDensitySaturation <: SaturationMethod
    ChemPotDensitySaturation(;vl = nothing,
                            vv = nothing,
                            crit = nothing,
                            f_limit = 0.0,
                            atol = 1e-8,
                            rtol = 1e-12,
                            max_iters = 10^4)
Saturation method for `saturation_pressure`. It uses equality of Chemical Potentials with a density basis. If no volumes are provided, it will use  [`x0_sat_pure`](@ref). 
`vl`  and `vl` are initial guesses for the liquid and vapour volumes.
`f_limit`, `atol`, `rtol`, `max_iters` are passed to the non linear system solver.
"""
function ChemPotDensitySaturation(;vl = nothing,
    vv = nothing,
    crit = nothing,
    f_limit = 0.0,
    atol = 1e-8,
    rtol = 1e-12,
    max_iters = 10^4)

    if (vl === nothing) && (vv === nothing)
        return ChemPotDensitySaturation{Nothing,typeof(crit)}(nothing,nothing,crit,f_limit,atol,rtol,max_iters)
    elseif !(vl === nothing) && (vv === nothing)
        vl = float(vl)
        return ChemPotDensitySaturation(vl,vv,crit,f_limit,atol,rtol,max_iters)
    elseif (vl === nothing) && !(vv === nothing)
        vv = float(vv)
        return ChemPotDensitySaturation(vl,vv,crit,f_limit,atol,rtol,max_iters)
    else
        T = one(vl)/one(vv)
        vl,vv,_ = promote(vl,vv,T)
        return ChemPotDensitySaturation(vl,vv,crit,f_limit,atol,rtol,max_iters)
    end
end

function saturation_pressure_impl(model::EoSModel, T, method::ChemPotDensitySaturation{Nothing})
    x0 = x0_sat_pure(model,T)
    vl,vv = x0
    method = ChemPotDensitySaturation(;vl,vv,method.crit)
    return saturation_pressure_impl(model,T,method)
end

function saturation_pressure_impl(model::EoSModel,T,method::ChemPotDensitySaturation)
    crit = method.crit
    if crit !== nothing
        Tc,_,_ = crit
        if Tc < T
            nan = zero(T)/zero(T)
            return (nan,nan,nan)
        end
    end
    return psat_chempot(model,T,method.vl,method.vv,NEqOptions(method))
end

function psat_chempot(model,T,vol_liq0,vol_vap0,options = NEqOptions())
    ρl0 = 1/vol_liq0
    ρv0 = 1/vol_vap0
    ρ0 = vec2(ρl0,ρv0,T)
    ofpsat = fobj_psat!(model, T)
    # sol = NLsolve.nlsolve(only_fj!(ofpsat), ρ0, method = :newton)
    sol = Solvers.nlsolve(ofpsat, ρ0, LineSearch(Newton()),options)
    #@show sol
    ρ = Solvers.x_sol(sol)
    vol_liq, vol_vap = 1 ./ ρ
    P = pressure(model, vol_vap, T)
    converged = check_valid_sat_pure(model,P,vol_liq,vol_vap,T)
    if converged
        return P, vol_liq, vol_vap
    else
        nan = zero(P)/zero(P)
        return (nan,nan,nan)
    end
end

export ChemPotDensitySaturation