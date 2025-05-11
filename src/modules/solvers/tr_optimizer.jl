include("tr_subproblem/TRS.jl")

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
            X,info = Solvers.TRS.trs_small(h, g, d, compute_local=true)
            s = X[1:n]
			outside = cnst(x, s)
			if d < sqrt(dmin)
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