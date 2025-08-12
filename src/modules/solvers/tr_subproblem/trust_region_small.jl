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

function check_interior!(X::AbstractMatrix, info::TRSinfo, P, q::AbstractVector; direct=false)
	if info.λ[1] < -1e-9 # Global solution is in the interior
		x1 = X[:, 1]
		if !direct
			cg!(x1, P, -q, tol=(eps(real(eltype(q)))/2)^(2/3))
		else
			x1 = -(P\q)
		end
		X[:, 1] = x1
		info.λ[1] = 0
	end
	if size(X, 2) > 1 && info.λ[2] < -1e-9
		# No local-no-global minimiser can exist in the interior; discard it.
		X = X[:, 1:1]; info.λ = info.λ[1:1]
	end
	return X, info
end
