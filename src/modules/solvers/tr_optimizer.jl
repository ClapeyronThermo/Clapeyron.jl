include("tr_subproblem/TRS.jl")

function Gershgorin(A)
    n = size(A)[1]
    e = Vector{eltype(A)}(undef,n)
    for j = 1:n
	sum = zero(eltype(A))
	for i = 1:n
	    sum += abs(A[i,j]);
	end
	e[j] = A[j,j] - (sum - abs(A[j,j]));
    end
    return minimum(e);
end

function exactstep(g, h, delta, tol, verbose)
    n = size(g)[1]
    max_iter = 10*n
    hard_case = false
	p = Vector{eltype(g)}(undef,n)
	q = Vector{eltype(g)}(undef,n)

    epsilon = eps(Float64)
    numeric_eps = abs(h[1,1])*epsilon
    for i = 2:n
    	if (abs(h[i,i])*epsilon > numeric_eps) 
    	    numeric_eps = abs(h[i,i])*epsilon
    	end
    end
    
    lambda = 0.0
    # We need to put a bracket around lambda.  It can't go below 0.  We
    # can get an upper bound using Gershgorin disks.  
    # This number is a lower bound on the eigenvalues in h
    
    h_min_eigenvalue = Gershgorin(h)

    normg = norm(g,2)
    lambda_min = 0.0
    lambda_max = normg / delta - h_min_eigenvalue
    hp = Matrix{Float64}(undef, n, n)
    hp = h
    
    # make hp symmetric if required maybe numerical unsymmetric
    symm = issymmetric(hp)

    # make hp symmetric
    if (!symm)
        for i = 2:n
            for j = 1:i-1
                hp[j,i] = h[i,j]
            end
	end
	symm = issymmetric(hp)

    end 
    iter = 0
    while (iter < max_iter)
        
        for i = 1:n
            hp[i,i] = h[i,i] + lambda
	end
	
	posdef = isposdef(hp)
	
	if (!posdef)
	    # If h is indefinite and g is equal to 0 then we should
	    # quit this loop and go right to the eigenvalue decomposition method.
	    if (normg <= numeric_eps)
	        hard_case = true
	        break
	    end

	    # narrow the bracket on lambda.  Obviously the current lambda is
	    # too small.

	    # jump towards the max value.  Eventually there will
	    # be a lambda that results in a cholesky decomposition.
	    alpha = 0.10
	    lambda += alpha*(lambda_max - lambda_min)
	    continue
	end

	# C = LU
    C = cholesky(hp)
    L = C.L
	p = g
	p = -p
    # Solve LU*p = -g for p.
	p = C \ p
	normp = norm(p,2)
	# Solve L*q = p for q.
	q = L \ p
	normq = norm(q,2)

	# check if we are done.  
	if (iter == 0 && posdef)
	    if (normp < delta)
	        # i will always be 0 in this case.
	        step_found = true
	        return p,h,iter,step_found,hard_case
	    end
	else
	    # if we are close enough to the solution then terminate
	    if (abs(normp - delta) < tol*delta)
		step_found = true
		h = hp
		return p,h,iter,step_found,hard_case
	    end
	end

	# shrink our bracket on lambda
	if (normp < delta)
	    lambda_max = lambda
	else
	    lambda_min = lambda
	end

	if (normp <= delta*epsilon)
	    alpha = 0.01
	    lambda = (1.0 - alpha)*lambda_min + alpha*lambda_max
	    continue
	end

    # figure out which lambda to try next
	lambda = lambda + normp*normp/normq/normq*(normp - delta) / delta

	# if we are outside the box we just bisect
	if (lambda <= lambda_min || lambda >= lambda_max)
	    lambda = 0.5*(lambda_min + lambda_max)
	end
	iter += 1
    end
    step_found = false
    return p,h,iter,step_found,hard_case
end

function trustregion_model(x, g, h)
    n = size(x)[1]
    hx = Vector{Float64}(undef,n)
    hx = h*x
    m = dot(x,g) + 0.5*dot(x,hx)
    return -m
end

function delta_Nocedal(func::Function, proj::Function, cnst::Function, d, dmin, dmax, x, lb, ub, f, g, h, verbose)
    c1 = 0.25
    c2 = 0.75
    c3 = 0.25
    c4 = 2.00
    eps = 2.2e-16
    
    n = size(x)[1]
    xp = Vector{Float64}(undef,n)
    
    # exactstep updates h to ensure its positive definition we use this updated hessian as our model
#    s,h,iter,step_found,hard_case = exactstep(g, h, d, 0.001,verbose)
    X,info = Solvers.TRS.trs_small(h, g, d, compute_local=true)
    s = X[1:n]
    
	use_backtracking = false
	if use_backtracking
		tau = 1.0
		outside = cnst(x, s)
		# our step is deemed to big to be feasible
		# we try to back track but if the step is to big we may not even take one as for example we are close to a bound
		# hence we try to tack a small step and the projection function will help us stay feasible
		while outside
			tau *= 0.5
			s .*= tau
			outside = cnst(x, s)
			if tau < sqrt(dmin)
				step_found = false
				break
			end
		end
	else
		outside = cnst(x, s)
		# our step is deemed to big to be feasible
		# we try to reduce trust region but if the step is to big we may not even take one as for example we are close to a bound
		# hence we try to tack a small step and the projection function will help us stay feasible
		while outside
		    d *= 0.5
			# exactstep updates h to ensure its positive definition we use this updated hessian as our model
		#    s,h,iter,step_found,hard_case = exactstep(g, h, d, 0.001, verbose)
            X,info = Solvers.TRS.trs_small(h, g, d, compute_local=true)
            s = X[1:n]
			outside = cnst(x, s)
			if d < 4.0*dmin
				step_found = false
				break
			end
		end
	end
    
    step = norm(s,2)
    
    xp = x .+ s
    xp = proj(xp)
    fp = func(xp)
    m = trustregion_model(s, g, h)
    rho = (f - fp) / abs(m)
    
    if rho < c1
        d = c3*d
    end
    
    # we use step >= d here as exactstep approaches d from above
    if rho > c2 && step >= d
        d = c4*d
    end
    
    # seems we may land on a saddle point accepting the point allows as to move away and start again
    # we accept the step even if rho <= 0.0. This causes d to approach dmin and exit if we are the best we can do
    
    x = xp
    f = fp
    
    d = max(d, dmin)
    d = min(d, dmax)
    
    return d,x,f

end

function trustregion_Dennis_Schnabel(	func::Function,
										grad::Function,
										hess::Function,
										proj::Function,
										cnst::Function,x,lb,ub,max_iters,tol,verbose)
    dmin = tol
    dmax = 0.5
    n = size(x)[1]
    check = false
    g = Vector{Float64}(undef,n)
    s = Vector{Float64}(undef,n)
    p = Vector{Float64}(undef,n)
    h = Matrix{Float64}(undef, n, n)
    
    x = proj(x)
    f = func(x)
    g = grad(x)
    h = hess(x)
    h = (h + h')/2.0
    
    d = dmax
    iter = 0
    p = x
    s = -g
 
    p = p .- g
    p = p .- x
    error = norm(p, Inf)

    while error > tol
        iter += 1
        d,x,f = delta_Nocedal(func::Function, proj::Function, cnst::Function, d, dmin, dmax, x, lb, ub, f, g, h,verbose)
        if d == dmin
			if error < sqrt(eps(Float64))
				# maybe best we can do
				break
			end
	    	check = true
	    	break
        end
        g = grad(x)
        h = hess(x)
        h = (h + h')/2.0
        p = x
        p = p .- g
        p = p .- x
        error = norm(p,Inf)
        if verbose
        	println("Iteration: $(iter) / $(max_iters) error = $(error)")
        end
        if iter >= max_iters
            check = true
            break
        end  
    end
    return x,f,iter,error,check
end
