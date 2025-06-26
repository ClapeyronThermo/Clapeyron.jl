mutable struct TRSinfo
	hard_case::Bool  # Flag indicating if we are in the hard case
	niter::Int # Number of iterations in eigs
	nmul::Int  # Number of multiplication with P in eigs
	λ::Vector  # Lagrange Multiplier(s)

	function TRSinfo(hard_case::Bool, niter::Int, nmul::Int, λ::Vector)
		new(hard_case, niter, nmul, λ)
	end
end

function trs_boundary(P, q::AbstractVector{T}, r::T, C::AbstractMatrix{T}; kwargs...) where {T}
	check_inputs(P, q, r, C)
	return trs_boundary((nev; kw...) -> gen_eigenproblem(P, q, r, C, nev; kw...),
		   (λ, V; kw...) -> pop_solution!(λ, V, P, q, r, C; kw...); kwargs...)
end

function trs_boundary(P, q::AbstractVector{T}, r::T; kwargs...) where {T}
	check_inputs(P, q, r)
	return trs_boundary((nev; kw...) -> eigenproblem(P, q, r, nev; kw...),
		   (λ, V; kw...) -> pop_solution!(λ, V, P, q, r, I; kw...); kwargs...)
end

function trs_boundary(solve_eigenproblem::Function, pop_solution!::Function;
	compute_local=false, tol_hard=2e-7, kwargs...)

	if compute_local
		nev=2  # We will need the two rightmost eigenvalues.
	else
		nev=1  # We will only need the rightmost eigenvalue
	end

	l, V, niter, nmult = solve_eigenproblem(nev; kwargs...)
	X, λ, hard_case = pop_solution!(l, V; tol_hard=tol_hard) # Pop global minimizer(s).
	if !compute_local || hard_case
		return X, TRSinfo(hard_case, niter, nmult, λ)
	else
		X_local, λ_local, _ = pop_solution!(l, V; tol_hard=tol_hard) # Pop local-no-global minimizer.
		return [X X_local], TRSinfo(hard_case, niter, nmult, [λ; λ_local])
	end
end

function check_inputs(P, q::AbstractVector{T}, r::T) where {T}
	@assert(issymmetric(P), "The cost matrix must be symmetric.")
	@assert(eltype(P) == T, "Inconsistent element types.")
	@assert(size(P, 1) == size(P, 2) == length(q), "Inconsistent matrix dimensions.")
end

function check_inputs(P, q::AbstractVector{T}, r::T, C::AbstractMatrix{T}) where {T}
	check_inputs(P, q, r)
	@assert(issymmetric(C), "The norm must be defined by a symmetric positive definite matrix.")
end

function pop_solution!(λ, V, P, q::AbstractVector{T}, r::T, C; tol_hard, direct=false) where {T}
	idx = argmax(real(λ)) 	# Pop rightmost eigenvector
	is_first_pop = all(isfinite.(λ))
	n = length(q)
	if abs(imag(λ[idx])) >= 1e-9 && !is_first_pop # A solution exists only on real eigenvalues
		return zeros(T, n, 0), zeros(T, 0), false
	end
	l = real(λ[idx]);
	complex_v = view(V, :, idx)
	if norm(real(complex_v)) > norm(imag(complex_v))
		v = real(complex_v)
	else
		v = imag(complex_v)  # Sometimes the retuned eigenvector is complex
	end
	v ./= norm(v)
	v1 = view(v, 1:n); v2 = view(v, n+1:2*n)

	if norm(q) < 1e-9
		# Solution the eigenvector of Px = lCx, which is parallel to v2
		x = r*v2/sqrt(v2'*C*v2)
		return [x -x], [l; l], true
	end

	# Extract solution
	norm_v1 = sqrt(dot(v1, C*v1))
	X = zeros(T, n, 0)
	hard_case = false
	if norm_v1 >= tol_hard
		if norm_v1 >= 1e-11
			x = -r*v1/norm_v1
			if sign(q'*v2) < 0
				x .= -x
			end

			residual = (P*x + l*C*x + q)
			# Polish the solution if necessary.
			if norm(residual) > 1e-7
				D = LinearMap{T}((x) -> P*x + l*C*x, n; issymmetric=true)
				lsqr!(x, D, -q, maxiter=20)
				x /= sqrt(x'*C*x)/r
			end
			X = reshape(x, n, 1)
		end
	elseif is_first_pop # hard case
		hard_case = true
		if !direct
			X = extract_solution_hard_case(P, q, r, C, l, v1, v2, tol_hard)
		else
			X = extract_solution_hard_case_direct(P, q, r, C, l, v1, v2)
		end
		# Sometimes extract_solution_hard_case fails in this case only return the Lagrange multipliers
		if !isreal(X)
			X = zeros(T, n, 0)
		end
	end
	λ[idx] = -Inf  # This ensures that the next pop_solution! would not consider the same solution.

	return X, l*ones(size(X, 2)), hard_case
end

function extract_solution_hard_case(P, q::AbstractVector{T}, r::T, C, l::T,
	v1::AbstractVector{T}, v2::AbstractVector{T}, tol::T) where T

	n = length(q)
	D = LinearMap{T}((x) -> P*x + l*C*x + v2*dot(v2, x)/100, n; issymmetric=true)
	y = minres(D, -q, maxiter=100, verbose=false)
	residual = P*y + l*C*y + q
	if sqrt(y'*C*y) > r || norm(residual) > tol
		l += T(1e-10) # Usually perturbing the multiplier and solving again works
		minres!(y, D, -q, tol=norm(residual)/tol, verbose=false)
	end

	if sqrt(y'*C*y) < r
		α = roots(Polynomial([y'*(C*y) - r^2, 2*(C*v2)'*y, v2'*(C*v2)]))
		x1 = y + α[1]*v2
		x2 = y + α[2]*v2
		return [x1 x2]
	else
		x = r*y/sqrt(y'*C*y)
		return reshape(x, n, 1)
	end
end

function extract_solution_hard_case_direct(P, q::AbstractVector{T}, r::T, C, l::T,
	v1::AbstractVector{T}, v2::AbstractVector{T}) where T

	n = length(q)
	λ, V = eigen(Symmetric(P + l*C))
	λ[abs.(λ) .< 1e-9] .= 0
	A = V*diagm(0 => λ)*V' # Essentially P + l*C with a "refined" nullspace
	F = qr(A, Val(true))
	y = -(F\q) # y is the minimum norm solution of P + l*C = q
	if sqrt(y'*C*y) < r
		α = roots(Polynomial([y'*(C*y) - r^2, 2*(C*v2)'*y, v2'*(C*v2)]))
		x1 = y + α[1]*v2
		x2 = y + α[2]*v2
		return [x1 x2]
	else
		x = r*y/sqrt(y'*C*y)
		return reshape(x, n, 1)
	end
end
