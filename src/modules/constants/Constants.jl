module Constants

using ForwardDiff

#gas constant
@inline R̄(::Type{ForwardDiff.Dual{V,T}}) where {V,T} = R̄(T)
R̄(::Type{Float64}) = 8.31446261815324
R̄(::Type{Float32}) = 8.314463f0
R̄(::Type{Float16}) = Float16(8.31f0)
R̄(::Type{BigFloat}) = big"8.31446261815324"
R̄(::Type{T}) where T = R̄(Float64)

#avogadro constant

@inline N_A(::Type{ForwardDiff.Dual{V,T}}) where {V,T} = N_A(T)
N_A(::Type{Float64}) = 6.02214076e23
N_A(::Type{Float32}) = 6.0221406f23
N_A(::Type{Float16}) = throw(DomainError("Avogadro constant cannot be represented in Float16 values"))
N_A(::Type{BigFloat}) = big"6.02214076e23"
N_A(::Type{T}) where T = N_A(Float64)

end #Constants module