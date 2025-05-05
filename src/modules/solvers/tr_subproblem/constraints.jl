function trs_boundary(P, q::AbstractVector{T}, r::T, A::AbstractMatrix{T}, b::AbstractVector{T}; kwargs...) where {T}
	project!, F = generate_nullspace_projector(A)
	x = find_feasible_point(b, r, project!, F)
	λ_max = max(maximum(eigs(P, nev=1, which=:LR)[1]), 0)
	return trs_boundary(P, q, r, project!, x, λ_max; kwargs...)
end

function generate_nullspace_projector(A::AbstractMatrix{T}) where T
	"""
	Generates a function project!(x) That projects x into the nullspace of A
	Also returns the factorization of the KKT matrix[I A'; A 0]
	"""
	F = factorize([I A'; A 0*I]) # KKT matrix
	n = size(A, 2)
	x_ = zeros(size(F, 1))
	y_ = zeros(size(F, 1))
	function project!(x::AbstractVector{T})
		copyto!(view(x_, 1:n), x)
		ldiv!(y_, F, x_)
		copyto!(x, view(y_, 1:n))
	end
	return project!, F
end

function find_feasible_point(b::AbstractVector{T}, r::T, project!, F::Factorization{T}) where T
	n = size(F, 1) - length(b)
	x = (F\[zeros(n); b])[1:n] # x is the minimizer of ‖x‖ with Ax = b
	@assert(norm(x) <= r, "The problem is infeasible.")
	d = project!(randn(n)) # Find a direction in the nullspace of A
	# Calculate alpha such that ‖x + alpha*d‖ = r
	alpha = roots(Polynomial([norm(x)^2 - r^2, 2*d'*x, norm(d)^2]))
	x += alpha[1]*d # Now ‖x‖ = r

	return x
end

function trs_boundary(P, q::AbstractVector{T}, r::T, project!, x::AbstractVector{T}, λ_max::T; kwargs...) where {T}
	"""
	Solves the TRS problem
	minimize    ½x'Px + q'x
	subject to  ‖x‖ = r
                Ax = b.
	Instead of passing A and b it is required to pass project! and x, where:
	- project!(x) projects (inplace) x to the nullspace of A; and
	- x is a point with ‖x‖ = r and Ax = b
	"""
	n = length(q)
	function p(y::AbstractVector{T}, x::AbstractVector{T}) where {T}
		mul!(y, P, x)
		# Substracting λ_max makes P negative definite which helps eigs for the constrained case
		axpy!(-λ_max, x, y)
		project!(y) # Project to the nullspace of A
	end
	x0 = x - project!(copy(x)) # x0 is perpendicular to the nullspace of A
	P_ = LinearMap{T}(p, n; ismutating=true, issymmetric=true)
	q_ = project!(q + P*x0 - λ_max*x0)
	r_ = norm(x - x0)
	output = trs_boundary((nev; kw...) -> eigenproblem(P_, q_, r_, nev;
			v0=[project!(randn(n)); project!(randn(n))], kw...,),
			(λ, V; kw...) -> pop_solution!(λ, V, P_, q_, r_, I; kw...); kwargs...)
	return shift_output(output..., x0, λ_max)
end

function shift_output(X, info, x0, λ_max)
	X .+= x0
	info.λ .-= λ_max
	return X, info
end

function trs(P, q::AbstractVector{T}, r::T, A::AbstractMatrix{T}, b::AbstractVector{T}; kwargs...) where {T}
	project!, F = generate_nullspace_projector(A)
	x = find_feasible_point(b, r, project!, F)
	λ_max = max(maximum(eigs(P, nev=1, which=:LR)[1]), 0)
	X, info = trs_boundary(P, q, r, project!, x, λ_max; kwargs...)

	# x0 is perpendicular to the nullspace of A
	x0 = X[:, 1] - project!(copy(X[:, 1]))
	# We remove x0 from the output(s) so that they belong to the nullspace of A
	shift_output(X, info, -x0, 0)
	# Define projected matrices (on the nullspace of A) for the cg
	P_projected = LinearMap{T}((y, x) -> project!(mul!(y, P, x)), length(x0);
		ismutating=true, issymmetric=true)
	q_projected = project!(q + P*x0)
	# Check interior solutions (if necessary) and add back x0
	X, info = check_interior!(X, info, P_projected, q_projected)
	return shift_output(X, info, x0, 0)
end