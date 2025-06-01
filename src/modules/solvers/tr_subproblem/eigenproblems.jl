function eigenproblem_small(P::AbstractMatrix, q::AbstractVector{T}, r::T, nev=1) where {T}
	"""
	Calculates eigenvalues/vectors of

	|-P    qq'/r^2|  |v1|  =  λ |v1|
	| I         -P|  |v2|       |v2|

	with eigen.
	"""
	n = length(q)
	A = zeros(T, 2*n, 2*n) # The matrix to perform eigendecomposition
	@inbounds for i in 1:n
		for j in 1:n
			A[i, j] = -P[i, j]
		end
		A[i + n, i] = one(T)
		for j in 1:n
			A[i + n, j + n] = -P[i, j]
		end
	end
	@inbounds for i in 1:n
		c = q[i]/r^2
		for j in n+1:2*n
			A[i, j] = c*q[j-n]
		end
	end
	λ, V = eigen!(A)
	return λ, V, 0, 0
end

function gen_eigenproblem_small(P::AbstractMatrix, q::AbstractVector{T}, r::T, C::AbstractMatrix, nev=1) where {T}
	"""
	Calculates rightmost eigenvalues/vectors of

	|-P    qq'/r^2|  |v1|  =  λ |C   0| |v1|
	| C         -P|  |v2|       |0   C| |v2|

	with Arpack's eigs.
	"""
	n = length(q)
	A = zeros(T, 2*n, 2*n) # Left matrix of eigenproblem
	@inbounds for i in 1:n
		for j in 1:n
			A[i, j] = -P[i, j]
		end
		for j in 1:n
			A[i + n, j] = C[i, j]
		end
		for j in 1:n
			A[i + n, j + n] = -P[i, j]
		end
	end
	@inbounds for i in 1:n
		c = q[i]/r^2
		for j in n+1:2*n
			A[i, j] = c*q[j-n]
		end
	end
	B = zeros(T, 2*n, 2*n) # Right matrix of eigenproblem
	@inbounds for i in 1:n
		for j = 1:n
			B[i, j] = C[i, j]
		end
		for j in 1:n
			B[i + n, j + n] = C[i, j]
		end
	end
	λ, V = eigen!(A, B)
	return λ, V, 0, 0
end
