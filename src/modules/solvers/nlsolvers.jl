struct BoundedLineSearch{B,LS} <: NLSolvers.LineSearcher
    lo::B
    hi::B
    ls::LS #actual line search
end

function find_alpha_bounds(x, d, lo, hi, alpha0 = one(Base.promote_eltype(x,d,lo,hi)))
    _1 = one(Base.promote_eltype(x,d,lo,hi))
    alpha = _1*alpha0
    _0 = zero(_1)
    for i in eachindex(x)
        if d[i] > _0
            # Moving toward upper bound
            alpha_i = 0.99*_1*(hi[i] - x[i]) / d[i]
            alpha = min(alpha, alpha_i)
        elseif d[i] < _0
            # Moving toward lower bound
            alpha_i = 0.99*_1*(lo[i] - x[i]) / d[i]
            alpha = min(alpha, alpha_i)
        end
        # If d[i] == 0, no constraint from this component
    end
    
    return alpha
end

BoundedLineSearch(lb,ub) = BoundedLineSearch(lb,ub,NLSolvers.Backtracking())

function NLSolvers.find_steplength(mstyle, ls::BoundedLineSearch{F,LS}, φ, λ::T) where {F,LS,T}
    d = φ.d
    x = φ.x
    λr = find_alpha_bounds(φ.x, φ.d, ls.lo, ls.hi, λ)
    if ls.ls isa NLSolvers.Static
        λx = T(λr*ls.ls.α)
        return λx, φ(λx), true
    end
    return NLSolvers.find_steplength(mstyle, ls.ls, φ, λr)
end

export BoundedLineSearch
