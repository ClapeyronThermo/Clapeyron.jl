function trs_small(P::AbstractMatrix{T}, q::AbstractVector{T}, r::T; kwargs...) where {T}
	output = trs_boundary_small(P, q, r; kwargs...)
	return check_interior!(output..., P, q; direct=true)
end

function trs_small(P::AbstractMatrix{T}, q::AbstractVector{T}, r::T, C::AbstractMatrix{T}; kwargs...) where {T}
	output = trs_boundary_small(P, q, r, C; kwargs...)
	return check_interior!(output..., P, q; direct=true)
end

function trs_boundary_small(P::AbstractMatrix{T}, q::AbstractVector{T}, r::T, C::AbstractMatrix{T}; kwargs...) where T
	check_inputs(P, q, r, C)
	return trs_boundary((nev; kw...) -> gen_eigenproblem_small(P, q, r, C, nev; kw...),
		   (λ, V; kw...) -> pop_solution!(λ, V, P, q, r, C; kw..., direct=true); kwargs...)
end

function trs_boundary_small(P::AbstractMatrix{T}, q::AbstractVector{T}, r::T; kwargs...) where {T}
	check_inputs(P, q, r)
	return trs_boundary((nev; kw...) -> eigenproblem_small(P, q, r, nev; kw...),
		   (λ, V; kw...) -> pop_solution!(λ, V, P, q, r, I; kw..., direct=true); kwargs...)
end