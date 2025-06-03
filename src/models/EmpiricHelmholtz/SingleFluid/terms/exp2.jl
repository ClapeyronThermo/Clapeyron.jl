mutable struct DoubleExpTerm  <: MultiParameterTerm
    n::Vector{Float64}
    t::Vector{Float64}
    d::Vector{Float64}
    gt::Vector{Float64}
    lt::Vector{Float64}
    gd::Vector{Float64}
    ld::Vector{Float64}
    function DoubleExpTerm(n,t,d,gt,lt,gd,ld)
        @assert length(n) == length(t) == length(d) == length(lt)
        @assert length(lt) == length(ld) == length(gt) == length(gd)
        return new(n,t,d,gt,lt,gd,ld)
        return param
    end
end

DoubleExpTerm() = DoubleExpTerm(Float64[],Float64[],Float64[],Float64[],Float64[],Float64[],Float64[])

active_term(term::DoubleExpTerm) = !iszero(length(term.n))

function Base.empty!(term::DoubleExpTerm)
    Base.empty!(term.n)
    Base.empty!(term.t)
    Base.empty!(term.d)
    Base.empty!(term.lt)
    Base.empty!(term.gt)
    Base.empty!(term.ld)
    Base.empty!(term.gd)
    return term
end

function a_term(term::DoubleExpTerm,δ,τ,lnδ,lnτ,_0)
    if active_term(term)
        n = term.n
        t = term.t
        d = term.d
        ld = term.ld
        gd = term.gd
        lt = term.lt
        gt = term.gt
        αᵣ = term_ar_exp2(δ,τ,lnδ,lnτ,_0,n,t,d,ld,gd,lt,gt)
    else
        αᵣ = _0
    end
    return αᵣ
end
