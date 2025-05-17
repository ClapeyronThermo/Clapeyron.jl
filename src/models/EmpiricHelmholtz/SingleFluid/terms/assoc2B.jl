mutable struct Associating2BTerm  <: MultiParameterTerm
    epsilonbar::Float64
    kappabar::Float64
    a::Float64
    m::Float64
    vbarn::Float64
end

Associating2BTerm() = Associating2BTerm(0.0,0.0,0.0,0.0,0.0)

active_term(term::Associating2BTerm) = !iszero(term.kappabar)

function Base.empty!(term::Associating2BTerm)
    term.epsilonbar = 0.0
    term.kappabar = 0.0
    term.a = 0.0
    term.m = 0.0
    term.vbarn = 0.0
    return term
end

function a_term(term::Associating2BTerm,δ,τ,lnδ,lnτ,_0)
    if active_term(term)
        ε̄ = term.epsilonbar
        κ̄ = term.kappabar
        a = term.a
        m = term.m
        v̄ₙ = term.v̄ₙ
        αᵣ = term_ar_assoc2b(δ,τ,lnδ,lnτ,_0,ε̄,κ̄,a,m,v̄ₙ)
    else
        αᵣ = _0
    end
    return αᵣ
end

    

