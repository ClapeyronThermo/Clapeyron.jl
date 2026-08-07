struct RestrictedLineSearch{F,LS} <: NLSolvers.LineSearcher
    f::F #function that restricts the line search
    ls::LS #actual line search
end

Base.@kwdef struct PositiveLinesearch{T,D,TOL}
    τ::T = 1.0
    decay::D =0.5
    tol::TOL = 1e-10
end

function NLSolvers.find_steplength(mstyle, ls::PositiveLinesearch, φ, λ::T) where T
    _d = φ.d
    _x = φ.z
    λx = T(positive_linesearch(_x, _d, λ; τ = ls.τ, decay = ls.decay, tol = ls.tol, s = -1))
    return λx, φ(λx), !isnan(λx)
end

function NLSolvers.find_steplength(mstyle, ls::RestrictedLineSearch{F,LS}, φ::T, λ) where {F,LS,T}
    λr = ls.f(φ,λ)
    NLSolvers.find_steplength(mstyle, ls.ls, φ, λ)
end


