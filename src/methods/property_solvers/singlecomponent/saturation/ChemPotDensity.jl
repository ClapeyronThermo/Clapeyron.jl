function ∂Helmholtz(model::EoSModel, ρ, T, z=[1.0])
    # Auxiliar functions to compute Helmholtz Energy and its first density derivative
    f(dρ) = eos(model, 1. / dρ, T, z)
    A, ∂A, = Solvers.f∂f(f, ρ)
    return A, ∂A
end

function ∂2Helmholtz(model::EoSModel, ρ, T, z=[1.0])
    # Auxiliar functions to compute Helmholtz Energy and its first and second density derivatives
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

    ps,μs =  scale_sat_pure(model)
    function f!(F,x)
        ρ_liq, ρ_vap = x
        A_liq, ∂A_liq = ∂Helmholtz(model, ρ_liq, T)
        A_vap, ∂A_vap = ∂Helmholtz(model, ρ_vap, T)
        P_liq = ∂A_liq*ρ_liq^2
        μ_liq = A_liq + ρ_liq * ∂A_liq
        P_vap = ∂A_vap*ρ_vap^2
        μ_vap = A_vap + ρ_vap * ∂A_vap
        F[1] = (μ_liq - μ_vap)μs
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


struct ChemPotDensitySaturation{T} <: SaturationMethod
    vl::Union{Nothing,T}
    vv::Union{Nothing,T}
end

"""
    ChemPotVSaturation()
    ChemPotDensitySaturation(;vl,vv)

Saturation method. It uses equality of Chemical Potentials with a density basis. If no volumes are provided, it will use  [`x0_sat_pure`](@ref). 

If those initial guesses fail and the specification is near critical point, it will try one more time, using Corresponding States instead.

`V0` is `[log10(Vₗ₀),log10(Vᵥ₀)]` , where `Vₗ₀`  and `Vᵥ₀` are initial guesses for the liquid and vapour volumes.
"""
function ChemPotDensitySaturation(;vl = nothing,vv = nothing)
    if (vl === nothing) && (vv === nothing)
        return ChemPotDensitySaturation{Nothing}(nothing,nothing)
    elseif !(vl === nothing) && (vv === nothing)
        vl = float(vl)
        return ChemPotDensitySaturation(vl,vv)
    elseif (vl === nothing) && !(vv === nothing)
        vv = float(vv)
        return ChemPotDensitySaturation(vl,vv)
    else
        T = one(vl)/one(vv)
        vl,vv,_ = promote(vl,vv,T)
        return ChemPotDensitySaturation(vl,vv)
    end
end

function saturation_pressure_impl(model::EoSModel, T, method::ChemPotDensitySaturation{Nothing})
    x0 = x0_sat_pure(model,T) .|> exp10
    vl,vv = x0
    method = ChemPotDensitySaturation(vl,vv)
    return saturation_pressure_impl(model,T,method)
end

function saturation_pressure_impl(model::EoSModel,T,method::ChemPotDensitySaturation)
    return psat_chempot(model,T,method.vl,method.vv)
end

function psat_chempot(model,T,vol_liq0,vol_vap0)
    ρl0 = 1/vol_liq0
    ρv0 = 1/vol_vap0
    ρ0 = vec2(ρl0,ρv0,T)
    ofpsat = fobj_psat!(model, T)
    # sol = NLsolve.nlsolve(only_fj!(ofpsat), ρ0, method = :newton)
    sol = Solvers.nlsolve(ofpsat, ρ0, LineSearch(Newton())) #LineSearch(Newton(),HZAW())
    #@show sol
    ρ = Solvers.x_sol(sol)
    vol_liq, vol_vap = 1 ./ ρ
    P = pressure(model, vol_vap, T)
    return P, vol_liq, vol_vap
end

export ChemPotDensitySaturation