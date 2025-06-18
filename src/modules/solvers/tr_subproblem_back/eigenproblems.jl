function eigenproblem(P, q::AbstractVector{T}, r::T, nev=1;
	tol=0.0, maxiter=300, v0=zeros((0,))) where {T}
	"""
	Calculates rightmost eigenvalues/vectors of

	|-P   qq'/r^2|  |v1|  =  λ|v1|
	| I        -P|  |v2|      |v2|

	with Arpack's eigs.
	"""
	n = length(q)
	if norm(q) < 1e-9
		# Solution is the minimum eigenvector of P
		(λ, V, nconv, niter, nmult, _) = eigs(P, nev=1, which=:SR, tol=tol, maxiter=maxiter)
		return -λ, [0*V; V], niter, nmult
	end

	function A(y::AbstractVector, x::AbstractVector)
		@inbounds y1 = view(y, 1:n); @inbounds y2 = view(y, n+1:2*n)
		@inbounds x1 = view(x, 1:n); @inbounds x2 = view(x, n+1:2*n)
		mul!(y1, P, x1)
		mul!(y2, P, x2)
		axpy!(-dot(q, x2)/r^2, q, y1)
		axpy!(-one(T), x1, y2)
	end
	D = LinearMap{T}(A, 2*n; ismutating=true)

	(λ, V, nconv, niter, nmult, _) = eigs(-D, nev=nev, which=:LR,
	tol=tol, maxiter=maxiter, v0=v0)
	@assert(nconv >= min(nev, 2), "Eigensolver has failed to converge.
		Try decreasing tolerance or changing parameters of eigs.")

	return λ, V, niter, 2*nmult
end

function gen_eigenproblem(P, q::AbstractVector{T}, r::T, C::AbstractMatrix, nev=1;
	tol=0.0, maxiter=500, v0=zeros((0,))) where {T}
	"""
	Calculates rightmost eigenvalues/vectors of

	|-P    qq'/r^2|  |v1|  =  λ |C   0| |v1|
	| C         -P|  |v2|       |0   C| |v2|

	with Arpack's eigs.
	"""
	n = length(q)
	function A(y::AbstractVector, x::AbstractVector)
		@inbounds y1 = view(y, 1:n); @inbounds y2 = view(y, n+1:2*n)
		@inbounds x1 = view(x, 1:n); @inbounds x2 = view(x, n+1:2*n)
		mul!(y1, C, x1)
		mul!(y2, P, x2)
		axpy!(-one(T), y1, y2)

		mul!(y1, P, x1)
		axpy!(-dot(q, x2)/r^2, q, y1)
	end
	D = LinearMap{T}(A, 2*n; ismutating=true)

	(λ, V, nconv, niter, nmult, _) = eigs(-D, [C 0*I; 0*I C], nev=nev, which=:LR, tol=tol)
	@assert(nconv >= min(nev, 2), "Eigensolver has failed to converge.
		Try decreasing tolerance or changing parameters of eigs.")

	return λ, V, niter, 2*nmult
end

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