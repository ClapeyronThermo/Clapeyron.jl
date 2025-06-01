function trs(P, q::AbstractVector{T}, r::T; kwargs...) where {T}
	output = trs_boundary(P, q, r; kwargs...)
	return check_interior!(output..., P, q)
end

function trs(P, q::AbstractVector{T}, r::T, C::AbstractMatrix{T}; kwargs...) where {T}
	output = trs_boundary(P, q, r, C; kwargs...)
	return check_interior!(output..., P, q)
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