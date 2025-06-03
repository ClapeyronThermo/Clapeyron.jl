struct NonAnalyticTerm  <: MultiParameterTerm
    A::Vector{Float64}
    B::Vector{Float64}
    C::Vector{Float64}
    D::Vector{Float64}
    a::Vector{Float64}
    b::Vector{Float64}
    beta::Vector{Float64}
    n::Vector{Float64}
    function NonAnalyticTerm(A,B,C,D,a,b,beta,n)
        @assert length(A) == length(B) == length(C) == length(D)
        @assert length(A) == length(a) == length(b) == length(beta)
        @assert length(beta) == length(n)
        return new(A,B,C,D,a,b,beta,n)
    end
end

NonAnalyticTerm() = NonAnalyticTerm(Float64[],Float64[],Float64[],Float64[],Float64[],Float64[],Float64[],Float64[])

active_term(term::NonAnalyticTerm) = !iszero(length(term.n))

function Base.empty!(term::NonAnalyticTerm)
    Base.empty!(term.A)
    Base.empty!(term.B)
    Base.empty!(term.C)
    Base.empty!(term.D)
    Base.empty!(term.a)
    Base.empty!(term.b)
    Base.empty!(term.beta)
    Base.empty!(term.n)
    return term
end

function a_term(term::NonAnalyticTerm,δ,τ,lnδ,lnτ,_0)
    if active_term(term)
        A,B,C,D,a,b,β,n = term.A,term.B,term.C,term.D,term.a,term.b,term.beta,term.n
        αᵣ = term_ar_na(δ,τ,lnδ,lnτ,_0,A,B,C,D,a,b,β,n)
    else
        αᵣ = _0
    end
    return αᵣ
end
