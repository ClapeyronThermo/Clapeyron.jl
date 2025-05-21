using JuMP, HiGHS

"""
    HELDTPFlash(;   max_HELD_iters::Int = 0
    				max_trust_region_iters::Int = 0
    				tol::Float64 = HELD_tol/10
    				HELD_tol::Float64 = sqrt(eps)
					add_pure_guess = true
    				add_anti_pure_guess = true
    				add_pure_component = [0]
    				add_random_guess = false
    				add_all_guess = false
    				verbose = true)

"""
Base.@kwdef struct HELDTPFlash <: TPFlashMethod
    max_HELD_iters::Int = 0
    max_trust_region_iters::Int = 0
    tol::Float64 = 0.1*sqrt(eps(Float64))
    HELD_tol::Float64 = sqrt(eps(Float64))
	add_pure_guess::Bool = true
    add_anti_pure_guess::Bool = true
    add_pure_component::Vector{Bool} = Vector{Bool}(undef,0)
    add_random_guess::Bool = false
    add_all_guess::Bool = false
    verbose::Bool = false
end

function tp_flash_impl(model::EoSModel, p, T, n, method::HELDTPFlash)
	z₀ = n
	sumz₀ = sum(z₀)
	z₀ ./= sumz₀
	if method.max_HELD_iters == 0
		max_HELD_iters = 100*length(n)
	else
		max_HELD_iters = method.max_HELD_iters
	end
	if method.max_trust_region_iters == 0
		max_trust_region_iters = 2000*length(n)
	else
		max_trust_region_iters = method.max_trust_region_iters
	end
	tol = method.tol
	HELD_tol = method.HELD_tol
	if length(method.add_pure_component) == 0 || length(method.add_pure_component) !== length(n)
		add_pure_component = fill(true,length(n))
	else
		add_pure_component = method.add_pure_component
	end
	add_pure_guess = method.add_pure_guess
	add_anti_pure_guess = method.add_anti_pure_guess
	add_random_guess = method.add_random_guess
	add_all_guess = method.add_all_guess
	verbose = method.verbose
	if verbose == true
		println("HELD  - Setup:")
		println("HELD  - trust region tolerence = $(tol)")
		println("HELD  - HELD tolerence = $(HELD_tol)")
		println("HELD  - add_pure_guess = $(add_pure_guess)")
		println("HELD  - add_anti_pure_guess = $(add_anti_pure_guess)")
		println("HELD  - add_pure_component = $(add_pure_component)")
		println("HELD  - add_random_guess = $(add_random_guess)")
		println("HELD  - add_all_guess = $(add_all_guess)")
	end
	beta,xp,vp,Gsol = HELD_impl(model,p,T,z₀,max_HELD_iters,max_trust_region_iters,tol,HELD_tol,add_pure_guess,add_anti_pure_guess,add_pure_component,add_random_guess,add_all_guess,verbose)

    return beta,xp,vp,Gsol
    
end

# new HELD

function Constraints(x,lb,ub,s)
    n = size(x)[1]
	xp = Vector{Base.promote_eltype(x,lb,ub,s)}(undef,n)
	
    outside = false

    for i=1:n
        xp[i] = x[i] + s[i]
    end
    for i = 1:n
       if xp[i] > ub[i]
	    	outside = true
	    end
		if xp[i] < lb[i]
	    	outside = true
		end
    end

    return outside

end

function ProjectionHELD(x,lb,ub)
    n = size(x)[1]
    p = Vector{Base.promote_eltype(x,lb,ub)}(undef,n)
    xp = Vector{Base.promote_eltype(x,lb,ub)}(undef,n)
    for i = 1:n-1
		p[i] = x[i]
		if p[i] < lb[i]
	    	p[i] = lb[i]
		end
		if p[i] > ub[i]
	    	p[i] = ub[i]
		end
		xp[i] = p[i]
    end
    sumxp = sum(xp[1:n-1])
    xp[n] = 1.0-sumxp
    if xp[n] < lb[1]
    	xp[n] = lb[1]
    end
    if xp[n] > ub[1]
    	xp[n] = ub[1]
    end
    sumxp = sum(xp[1:n])
    xp ./= sumxp
    for i = 1:n-1
    	p[i] = xp[i]
    end
    p[n] = x[n]
    if p[n] < lb[n]
		p[n] = lb[n]
	end
	if p[n] > ub[n]
	    p[n] = ub[n]
	end
    return p
end

function Projection(x,lb,ub)
    n = size(x)[1]
    p = Vector{Base.promote_eltype(x,lb,ub)}(undef,n)
    for i = 1:n
		p[i] = x[i]
		if p[i] < lb[i]
	    	p[i] = lb[i]
		end
		if p[i] > ub[i]
	    	p[i] = ub[i]
		end 
    end
    return p
end

function HELD_func_rho(model,p,T,x₀,vref,rho)
    v = vref/rho
    A = eos(model,v,T,x₀)
    f = (A + p*v)/R̄/T
    return f
end

function HELD_volume1(model,p,T,x₀,vref,rho)
	G(x) = HELD_func_rho(model,p,T,[x,1.0-x],vref,rho) + 5.474135924873183*(0.95 - x)
    dG(x) = Solvers.derivative(G,x)
	ddG(x) = Solvers.derivative(dG,x)
    return G(x₀), dG(x₀), ddG(x₀)
end

function HELD_volume2(model,p,T,x₀,vref,rho)
	G(x) = HELD_func_rho(model,p,T,x₀,vref,x)
    dG(x) = Solvers.derivative(G,x)
	ddG(x) = Solvers.derivative(dG,x)
    return G(rho), dG(rho), ddG(rho)
end

function HELD_volume(model,p,T,x₀)
	pure = split_pure_model(model)
	crit = crit_pure.(pure)
	vref = 0.0
	for i = 1:length(x₀)
    	Tc,pc,vc = crit[i]
    	vref += x₀[i]*vc
	end
	v₀ = HELD_density(model,p,T,x₀,vref)
	return vref/v₀
end

function HELD_density(model,p,T,x₀,vref)
	G(x) = HELD_func_rho(model,p,T,x₀,vref,x)
	dG(x) = Solvers.derivative(G,x)
	ddG(x) = Solvers.derivative(dG,x)
	
	# calculate rho_ideal
	rho_ideal = vref/(R̄*T/p)
	rho_min = 1.0e-6
	rho_max = 15.0

	drho = rho_ideal/2.0
	if drho < rho_min
		rho_min = drho
		drho = 2.0*rho_min
	else
		drho = (rho_min + rho_ideal)/2.0
	end

	rho1 = rho_min
	rho2 = rho1+drho
	if abs(dG(rho2)) < sqrt(eps(Float64))
		rho2 += drho
	end

	rho_bracket = Vector{Vector{Float64}}(undef,0)
	while rho2 <= rho_max
		if ddG(rho1)*ddG(rho2) < 0.0
			push!(rho_bracket,[rho1,rho2])
		end
		if rho1 >= 0.01
			drho = 0.01
		end
		rho1=rho2
		rho2=rho1+drho
		if abs(ddG(rho2)) < sqrt(eps(Float64))
			rho2 += drho
		end
	end

#	println("rho_bracket  = $(rho_bracket)")

	rho_spinodial = Vector{Float64}(undef,0)
	for ib = eachindex(rho_bracket)
		ans = Roots.find_zero(ddG, rho_bracket[ib])
		for ia = eachindex(ans)
			push!(rho_spinodial,ans[ia])
		end
	end

	rho_sp_low = rho_min
	rho_sp_high = rho_max

	if length(rho_spinodial) > 1
		rho_sp_low = rho_spinodial[1]
		rho_sp_high = rho_spinodial[end]
	end

#	println("rho_spinodial  = $(rho_spinodial)")
#	println("rho_sp_low  = $(rho_sp_low)")
#	println("rho_sp_high  = $(rho_sp_high)")

	drho = rho_ideal/2.0
	if drho < rho_min
		rho_min = drho
		drho = 2.0*rho_min
	else
		drho = (rho_min + rho_ideal)/2.0
	end

	rho1 = rho_min
	rho2 = rho1+drho
	if abs(dG(rho2)) < sqrt(eps(Float64))
		rho2 += drho
	end

	rho_bracket = Vector{Vector{Float64}}(undef,0)
	while rho2 <= rho_max
		if dG(rho1)*dG(rho2) < 0.0
			push!(rho_bracket,[rho1,rho2])
		end
		if rho1 >= 0.01
			drho = 0.01
		end
		rho1=rho2
		rho2=rho1+drho
		if abs(dG(rho2)) < sqrt(eps(Float64))
			rho2 += drho
		end
	end

#	println("rho_bracket  = $(rho_bracket)")

	rho_found = Vector{Float64}(undef,0)
	for ib = eachindex(rho_bracket)
		ans = Roots.find_zero(dG, rho_bracket[ib])
		for ia = eachindex(ans)
			stab = ddG(ans[ia])
			if stab > 0.0
				push!(rho_found,ans[ia])
			end
		end
	end

#	println("rho_found  = $(rho_found)")

	rho_stable_set = Vector{Float64}(undef,0)
	if length(rho_spinodial) > 1
		if rho_found[1] < rho_sp_low
			push!(rho_stable_set,rho_found[1])
		end
		if rho_found[end] > rho_sp_high
			push!(rho_stable_set,rho_found[end])
		end
	else
		push!(rho_stable_set,rho_found[1])
		push!(rho_stable_set,rho_found[end])
	end

#	println("rho_stable_set  = $(rho_stable_set)")

	rho_stable = rho_stable_set[1]
	if length(rho_stable_set) > 1
		if G(rho_stable_set[end]) < G(rho_stable)
			rho_stable = rho_stable_set[end]
		end
	end

	return rho_stable

end

function HELD_func(model,p,T,n₀,v₀,x,λ)
    nc = length(n₀)
    xₙ = append!(deepcopy(x[1:nc-1]),1.0 - sum(x[1:nc-1]))
    v = v₀/x[end]
    A = eos(model,v,T,xₙ)
    f = (A + p*v)/R̄/T + ∑(λ.*(n₀[1:nc-1] .- xₙ[1:nc-1]))
    return f
end

function Gibbs_func(model,p,T,n₀,v₀, np, x)
	nc = length(n₀)
	beta = Vector{eltype(x)}(undef,0)
	xp = Vector{Vector{eltype(x)}}(undef,0)
    vp = Vector{eltype(x)}(undef,0)
    ix = 0
    for ip=1:np-1
    	ix += 1
    	push!(beta,deepcopy(x[ix]))
    	xc = Vector{eltype(x)}(undef,0)
    	for ic = 1:nc-1
    		ix += 1
    		push!(xc,deepcopy(x[ix]))
    	end
    	push!(xc,1.0 - sum(xc[1:nc-1]))
    	push!(xp,xc)
    	ix += 1
    	push!(vp,deepcopy(x[ix]))
    end
    ix += 1
    push!(vp,deepcopy(x[ix]))
    push!(beta,1.0 - sum(beta[1:np-1]))
    xc = Vector{eltype(x)}(undef,0)
    for ic = 1:nc
    	sum_moles = 0.0
    	for ip = 1:np-1
    		sum_moles += beta[ip]*xp[ip][ic]
    	end
    	push!(xc,(n₀[ic] - sum_moles)/beta[np])
    end
	push!(xp,xc)
    f = 0.0
    for ip=1:np
    	v = v₀/vp[ip]
	   	f += beta[ip]*(eos(model,v,T,xp[ip]) + p*v)/R̄/T
    end
    return f
end

function initial_compositions(model,p,T,z,add_pure_guess,add_anti_pure_guess,add_pure_component,add_random_guess,add_all_guess)
    
	n = length(z)

    lb = zeros(n)
    ub = ones(n)
    for i=1:n
        lb[i] += eps(Float64)*100.0
        ub[i] -= eps(Float64)*100.0
    end
	
    xp = Vector{Vector{Float64}}(undef,0)
           
    # Wilson k-values
    K = wilson_k_values(model,p,T)
    # vapour liquid like estimates
    xvap = fill(1.,n)
    xvap = z.*K
    sumxvap = sum(xvap)
    xvap ./= sumxvap
    xvap = Projection(xvap,lb,ub)
    sumxvap = sum(xvap)
    xvap ./= sumxvap
    push!(xp,xvap)
    xliq = fill(1.,n)
    xliq = z./K
    sumxliq = sum(xliq)
    xliq ./= sumxliq
    xliq = Projection(xliq,lb,ub)
    sumxliq = sum(xliq)
    xliq ./= sumxliq
    push!(xp,xliq)

   # Wilson k-values to power 1/3 for close to critical point
    Kn = K.^(1/3)
    xvap = z.*Kn
    sumxvap = sum(xvap)
    xvap ./= sumxvap
    xvap = Projection(xvap,lb,ub)
    sumxvap = sum(xvap)
    xvap ./= sumxvap
    push!(xp,xvap)
    xliq = z./Kn
    sumxliq = sum(xliq)
    xliq ./= sumxliq
    xliq = Projection(xliq,lb,ub)
    sumxliq = sum(xliq)
    xliq ./= sumxliq
    push!(xp,xliq)

   # Wilson k-values to power 3 for leaner vapour and richer liquid
    Kn = K.^3
    xvap = z.*Kn
    sumxvap = sum(xvap)
    xvap ./= sumxvap
    xvap = Projection(xvap,lb,ub)
    sumxvap = sum(xvap)
    xvap ./= sumxvap
    push!(xp,xvap)
    xliq = z./Kn
    sumxliq = sum(xliq)
    xliq ./= sumxliq
    xliq = Projection(xliq,lb,ub)
    sumxliq = sum(xliq)
    xliq ./= sumxliq
    push!(xp,xliq)
    
    if add_pure_guess || add_all_guess
  		# pure generation
    	k = 1000.0 - (n-1)
    	for i in 1:n
    		if add_pure_component[i] == true
        		xi = fill(1.,n)
        		xi[i] = k
        		sumxi = sum(xi)
        		xi ./= sumxi
        		push!(xp,xi)
        	end
    	end
    end
    
    if add_anti_pure_guess || add_all_guess
    	# anti pure generation
    	k = 1000.0 - (n-1)
    	for i in 1:n
    		if add_pure_component[i] == true
        		xi = z
        		xi[i] = z[i]/k
        		sumxi = sum(xi)
        		xi ./= sumxi
        		xi = Projection(xi,lb,ub)
    			sumxi = sum(xi)
    			xi ./= sumxi
        		push!(xp,xi)
        		# anti pure Wilson
        		K = wilson_k_values(model,p,T)
    			# vapour liquid like estimates
    			xvap = fill(1.,n)
    			xvap = xi.*K
    			sumxvap = sum(xvap)
    			xvap ./= sumxvap
    			xvap = Projection(xvap,lb,ub)
    			sumxvap = sum(xvap)
    			xvap ./= sumxvap
    			push!(xp,xvap)
    			xliq = fill(1.,n)
    			xliq = xi./K
    			sumxliq = sum(xliq)
    			xliq ./= sumxliq
    			xliq = Projection(xliq,lb,ub)
    			sumxliq = sum(xliq)
    			xliq ./= sumxliq
    			push!(xp,xliq)
    		end
    	end
    end
    
    if add_random_guess || add_all_guess
    	# some random guesses
    	for ir = 1:n
        	xr = fill(0.,n)
    		for i = 1:n
    			xr[i] = 1.0e-6 + (1.0 - 2e-6)*rand()
    		end
    		sumxr = sum(xr)
    		xr ./= sumxr
    		push!(xp,xr)
   		end
    end
    
    return xp
    
end

function Pereira_compositions(model,p,T,z)
    
	n = length(z)

    lb = zeros(n)
    ub = ones(n)
    for i=1:n
        lb[i] += eps(Float64)*100.0
        ub[i] -= eps(Float64)*100.0
    end
	
    xp = Vector{Vector{Float64}}(undef,0)
	d = Vector{Float64}(undef,0)

	use_full_dset = false
	if use_full_dset
		nd = 2.0^(n-1)
		Dd = 1.0/nd
		dd = 0.0
		for id = 1:nd-1
			dd += Dd
			push!(d,dd)
		end
	else
		push!(d,0.5)
	end
    
	for id in eachindex(d)
		for i = 1:n-1
			x̂ = fill(0.0,n)
			x̄ = fill(0.0,n)
			for j = 1:n
				if (j == i)
					x̂[j] = d[id]*z[i]
					x̄[j] = z[i] + d[id]*(1.0 - z[i])
				else
					x̂[j] = (1.0 - d[id]*z[i])/(n-1)
					x̄[j] = (1.0 - (z[i] + d[id]*(1.0 - z[i])))/(n-1)
				end
			end
			x̂ = Projection(x̂,lb,ub)
			sumx̂ = sum(x̂)
			x̂ ./= sumx̂
			push!(xp,x̂)
			x̄ = Projection(x̄,lb,ub)
			sumx̄ = sum(x̄)
			x̄ ./= sumx̄
			push!(xp,x̄)
		end
	end

    return xp
    
end

function HELD_impl(model,p,T,z₀,
	max_HELD_iters,
	max_trust_region_iters,
	tol,
	HELD_tol,
	add_pure_guess,
	add_anti_pure_guess,
	add_pure_component,
	add_random_guess,
	add_all_guess,
	verbose)
	
	# z₀ must sum to one, i.e. it is a mole fraction vector
    nc = length(z₀)
  # calculate reference volume based on Kays rule and vc[i] and scale to give water a vref/v ~ 1
    pure = split_pure_model(model)
    crit = crit_pure.(pure)
    vref = 0.0
    for i= 1:nc
    	Tc,pc,vc = crit[i]
    	vref += z₀[i]*vc
    end
 	ρ₀ = HELD_density(model,p,T,z₀,vref)
	v₀ = vref/ρ₀
    μ₀ = VT_chemical_potential(model,v₀,T,z₀)
    λ₀ = (μ₀[1:nc-1] .- μ₀[nc])/R̄/T
 
	ρ₀ = HELD_density(model,p,T,z₀,vref)
    G(x) = HELD_func(model,p,T,z₀,vref,x,λ₀)
    G_g(x) = Solvers.gradient(G,x)
    G_h(x) = Solvers.hessian(G,x)
    lb = zeros(nc)
    ub = ones(nc)
    for i=1:nc
        lb[i] += eps(Float64)*100.0
        ub[i] -= eps(Float64)*100.0
    end
    ub[nc] = 15.0
    projHELD(x) = ProjectionHELD(x,lb,ub)
	cnstHELD(x,s) = Constraints(x,lb,ub,s)
 	x₀ = append!(deepcopy(z₀[1:nc-1]),ρ₀)
    G₀ = G(x₀)
    if verbose == true
    		println("HELD Step 1 - Initialisation:")
    		println("HELD Step 1 - UBDⱽ = $(G₀)")
    		println("HELD Step 1 - λ₀ = $(λ₀)")
    end
    xi = initial_compositions(model,p,T,z₀,add_pure_guess,add_anti_pure_guess,add_pure_component,add_random_guess,add_all_guess)
    fmins = Vector{Float64}(undef,0)
    xmins = Vector{Vector{Float64}}(undef,0)
    for ix = 1:length(xi)
		ρi = HELD_density(model,p,T,xi[ix],vref)
		xρi = append!(deepcopy(xi[ix][1:nc-1]),ρi)
    	xmin,fmin,iter,error,check = Solvers.trustregion_Dennis_Schnabel(G, G_g, G_h, projHELD,cnstHELD, xρi, lb, ub, max_trust_region_iters, tol, false)
    
	#	if verbose == true
    #    	println("HELD Step 3 - IPₓᵥ solve, fmin = $(fmin) error = $(error) iter = $(iter)")
	#		println("HELD Step 3 - IPₓᵥ solve, xmin = $(xmin) check = $(check)")
    #	end 
    	
		# we add fmin < G₀ as we are searching for instability
    	if fmin < G₀ && check == false
    		push!(fmins,fmin)
    		push!(xmins,xmin)
    	end
    end
    
    fmins_unique, xmins_unique, stable = HELD_clean_local_solutions(G₀, x₀, fmins, xmins, tol, verbose)
    
    if verbose == true
    		println("HELD Step 1 - Phase stability check completed: $(length(fmins_unique)) unique solutions found")
    end
    if stable
        # return starting solution as it is stable
        if verbose == true
        	println("HELD Step 1 - Fluid is stable")
    	end
    	if verbose == true
			println("HELD Step 6 - Complete")
			println("HELD Step 6 - Phase found 1")
			println("HELD Step 6 - Phase moles:")
			println("HELD Step 6 - Phase beta(1) = 1")
			println("HELD Step 6 - Phase mole fraction:")
			println("HELD Step 6 - Phase x(1) = $(z₀)")
			println("HELD Step 6 - Phase volumes:")
			println("HELD Step 6 - Phase volume(1) = $(v₀)")
			println("HELD Step 6 - Minimum Gibbs Energy = $(G₀)")
    	end
    	return  [1.0], [z₀], [v₀], G₀
    else
        # return unique solutions, these need to be added to M and are good initial guesses for phases
        if verbose == true
            println("HELD Step 1 - Fluid is unstable, search for phases begins")
            println("HELD Step 1 - Initialise set ℳ used for OPₓᵥ , and ℳguess used for local minimisations in IPₓᵥ")
    	end
    	
        # for cutting plane we are working on the Gibbs surface so λ is zero
        λi = fill(0.,nc-1)
        Gi(x) = HELD_func(model,p,T,z₀,vref,x,λi)
        
        # set up inital ℳ used for OPₓᵥ
        # ℳi is [xi[1:nc-1], Vref/Vi, Gi]
    	ℳ = Vector{Vector{Float64}}(undef,0)
		xm = Pereira_compositions(model,p,T,z₀)
		# set up initial ℳ set
		# add initial guesses and the newly found minimums from first iteration stability check
		for im = 1:length(xm)
			ρm = HELD_density(model,p,T,xm[im],vref)
			xρm = append!(deepcopy(xm[im][1:nc-1]),ρm)
    		xρGim = append!(deepcopy(xρm),Gi(xρm))
    		push!(ℳ,xρGim)
    	end
		for ii = 1:length(xi)
			ρi = HELD_density(model,p,T,xi[ii],vref)
			xρi = append!(deepcopy(xi[ii][1:nc-1]),ρi)
    		xρGii = append!(deepcopy(xρi),Gi(xρi))
    		push!(ℳ,xρGii)
    	end
		for i = 1:length(fmins_unique)
			xminsGi_unique = append!(deepcopy(xmins_unique[i]),Gi(xmins_unique[i]))
			push!(ℳ,xminsGi_unique)
		end

        # set up inital ℳguess used for local minimisations in IPₓᵥ
        # ℳguessi is [xi[1:nc-1], Vref/Vi, Gi]
		ℳguess = Vector{Vector{Float64}}(undef,0)
		for ii = 1:length(xi)
			ρi = HELD_density(model,p,T,xi[ii],vref)
			xρi = append!(deepcopy(xi[ii][1:nc-1]),ρi)
			xρGii = append!(deepcopy(xρi),Gi(xρi))
    		push!(ℳguess,xρGii)
    	end
		for i = 1:length(fmins_unique)
			xminsGi_unique = append!(deepcopy(xmins_unique[i]),Gi(xmins_unique[i]))
			push!(ℳguess,xminsGi_unique)
		end
		
		n_unique_previous = length(fmins_unique)

		# now we have our first ℳ we can solve the cutting plane problem to get new λ
		# with new λ we can solve the local minimisation to get new unique minimums to add to ℳ
    	
    	UBDⱽ =  G₀
    	LBDⱽ = -Inf

		# need some way of getting the estimated lower and upper bounds on λ, 10 seems OK so far but may not be universal
		λnorm = norm(λ₀,Inf)
		λmax = 1.2*λnorm
		λᴸ = fill(-λmax,nc-1)
		λᵁ = fill( λmax,nc-1)

		limit_λs_by_bounds = true
    	
    	use_global_solution = false

    	λStalling_count = 0

    	np = 1
    	HELD_complete = false
    	xHELD = Vector{Float64}(undef,0)
    	
    	for k = 1:max_HELD_iters

        	if verbose == true
        		println("HELD Step 2 - iteration $(k)")
            	println("HELD Step 2 - Solve OPₓᵥ for new λˢ")
    		end
    		
    		OPₓᵥ = Model(HiGHS.Optimizer)
    		set_optimizer_attribute(OPₓᵥ, "log_to_console", false)
    		set_optimizer_attribute(OPₓᵥ, "output_flag", false)
    
    		@variable(OPₓᵥ, v)
    		@variable(OPₓᵥ, λ[1:nc-1])
    
    		@constraint(OPₓᵥ,v <= UBDⱽ)
    		#v <= Gi + ∑(λi*(n-xi))
    		@constraint(OPₓᵥ,[i ∈ 1:length(ℳ)],v <= ℳ[i][nc+1]+sum(λ.*(z₀[1:nc-1] .- ℳ[i][1:nc-1])))

			if limit_λs_by_bounds
    			#λᴸ <= λ <= λᵁ 
    			@constraint(OPₓᵥ,[i ∈ 1:nc-1],λᴸ[i] <= λ[i] <= λᵁ[i])
			end

    		@objective(OPₓᵥ, Max, v)
			use_ipm = false
			if use_ipm
    			ipm_tol = max(tol, 1.0e-10)
    			set_attribute(OPₓᵥ, "solver", "ipm")    		
    			set_attribute(OPₓᵥ, "ipm_optimality_tolerance", ipm_tol)
			else
				simplex_tol = max(tol, 1.0e-10)
	   			set_attribute(OPₓᵥ, "solver", "simplex")
				set_attribute(OPₓᵥ, "dual_feasibility_tolerance", simplex_tol)
				set_attribute(OPₓᵥ, "primal_feasibility_tolerance", simplex_tol)
				set_attribute(OPₓᵥ, "dual_residual_tolerance", simplex_tol)
				set_attribute(OPₓᵥ, "primal_residual_tolerance", simplex_tol)
			end
    		optimize!(OPₓᵥ)
    		λˢ = JuMP.value.(λ)

    		UBDⱽ  = JuMP.value.(v)

			λnorm = norm(λˢ,Inf)
			if k>nc
				filter = 0.1
				λmax  = filter*1.05*λnorm + (1.0 - filter)*λmax
			end
			λᴸ = fill(-λmax,nc-1)
			λᵁ = fill( λmax,nc-1)

			if verbose == true
    			println("HELD Step 2 - Update UBDⱽ and λˢ from OPₓᵥ: UBDⱽ = $(UBDⱽ)")
        		println("HELD Step 2 - λˢ = $(λˢ)")
				println("HELD Step 2 - λnorm = $(λnorm) λmax = $(λmax)")
    		end
    		
    		Dλ = λˢ  .- λ₀
			λStalling = false
    		if norm(Dλ,Inf) < HELD_tol
    		    if verbose == true
        			println("HELD Step 2 - λˢ has stalled use random inital guesses to provide a chance to converge")
    			end
				λStalling = true
				λStalling_count += 1
    		else
    			λStalling = false
    			λ₀ = λˢ 
    		end

    	    if verbose == true
        		println("HELD Step 3 - IPₓᵥ solve, generate cutting plane with λˢ")
    		end
    		Gˢ(x) = HELD_func(model,p,T,z₀,vref,x,λˢ)
    		Gˢ_g(x) = Solvers.gradient(Gˢ, x)
    		Gˢ_h(x) = Solvers.hessian(Gˢ, x)
    		
    		fmins = Vector{Float64}(undef,0)
    		xmins = Vector{Vector{Float64}}(undef,0)
    			
    		for ix = 1:length(ℳguess)
    			xmin,fmin,iter,error,check = Solvers.trustregion_Dennis_Schnabel(Gˢ, Gˢ_g, Gˢ_h, projHELD, cnstHELD, ℳguess[ix][1:nc], lb, ub, max_trust_region_iters, tol, false)
    			
   			#	if verbose == true
        	#		println("HELD Step 3 - IPₓᵥ solve, fmin = $(fmin) error = $(error) iter = $(iter)")
			#		println("HELD Step 3 - IPₓᵥ solve, xmin = $(xmin)")
 	  		#	end
				
				if fmin <= UBDⱽ  && check == false
					push!(fmins,fmin)
					push!(xmins,xmin)
				end

    		end
		
			fmins_unique, xmins_unique, stable = HELD_clean_local_solutions(UBDⱽ, x₀, fmins, xmins, tol, verbose)
			
			if length(fmins_unique) < 1 || λStalling
				use_global_solution = true
			else
				use_global_solution = false
			end
			
			if verbose == true
        		println("HELD Step 3 - IPₓᵥ solve, $(length(fmins_unique)) unique solutions found")
#        		println("HELD Step 3 - IPₓᵥ solve, xmins_unique = $(xmins_unique) unique solutions found")
    		end
    		
    		ℒ = Vector{Vector{Float64}}(undef,0)
    		if length(fmins_unique) > 0
				# find lowest minimum of returned set.
				LBDⱽ = fmins_unique[1]
				iLBDⱽ = 1
				for i = 2:length(fmins_unique)
					if fmins_unique[i] < LBDⱽ
				    	LBDⱽ = fmins_unique[i]
				    	iLBDⱽ = i
					end 
				end
				
				# sometimes the are more than one solution that can be added
				for i = 1:length(fmins_unique)
					if abs(fmins_unique[i] - fmins_unique[iLBDⱽ]) < tol
						push!(ℒ,xmins_unique[i])
					end
				end
				
				if verbose == true
        			println("HELD Step 3 - ℒ set contains $(length(ℒ)) items")
  				end
				
			end
			
			solution_found = true
			if use_global_solution
    			if verbose == true
        			println("HELD Step 3 - IPₓᵥ Global solution required")
  				end
				iter_rand_max = 10*nc*nc
				for iter_rand = 1:iter_rand_max
					xr = fill(0.,nc)
    				for i = 1:nc
    						xr[i] = 1.0e-6 + (1.0 - 2e-6)*rand()
    				end
    				sumxr = sum(xr)
    				xr ./= sumxr
    				xr = Projection(xr,lb,ub)
    				sumxr = sum(xr)
    				xr ./= sumxr
					ρr = HELD_density(model,p,T,xr,vref)
					xρr = append!(deepcopy(xr[1:nc-1]),ρr)
					xρGr = append!(deepcopy(xρr),Gi(xρr))
    				xmin,fmin,iter,error,check = Solvers.trustregion_Dennis_Schnabel(Gˢ, Gˢ_g, Gˢ_h, projHELD, cnstHELD,  xρGr[1:nc], lb, ub, max_trust_region_iters, tol, false)
    				
#    				if verbose == true
#        				println("HELD Step 3 - IPₓᵥ solve, fmin = $(fmin) error = $(error) iter = $(iter)")
#    				end

					if fmin <= UBDⱽ  && check == false
    					push!(fmins,fmin)
    					push!(xmins,xmin)
    				end
				end
				
				fmins_unique, xmins_unique, stable = HELD_clean_local_solutions(UBDⱽ, x₀, fmins, xmins, tol, verbose)
				
				if length(fmins_unique) > 0
					# find lowest minimum of returned set.
					LBDⱽ = fmins_unique[1]
					iLBDⱽ = 1
					for i = 2:length(fmins_unique)
						if fmins_unique[i] < LBDⱽ
				    		LBDⱽ = fmins_unique[i]
				    		iLBDⱽ = i
						end 
					end
					
					# sometimes the are more than one solution that can be added
					for i = 1:length(fmins_unique)
						if abs(fmins_unique[i] - fmins_unique[iLBDⱽ]) < tol
							push!(ℒ,xmins_unique[i])
						end
					end
				
					if verbose == true
        				println("HELD Step 3 - Global solution - ℒ set contains $(length(ℒ)) items")
  					end
  					
				else
					solution_found = false
				end	
    		end
    		
    		if !solution_found || λStalling_count > 3
    			# return starting solution as we have no phases to add to the solution so this is the best we can do
        		if verbose == true
        			println("HELD Step 1 - Global solution failed, its wise to check this solution")
    			end
    			if verbose == true
					println("HELD Step 6 - Complete")
					println("HELD Step 6 - Phase found 1")
					println("HELD Step 6 - Phase moles:")
					println("HELD Step 6 - Phase beta(1) = 1")
					println("HELD Step 6 - Phase mole fraction:")
					println("HELD Step 6 - Phase x(1) = $(z₀)")
					println("HELD Step 6 - Phase volumes:")
					println("HELD Step 6 - Phase volume(1) = $(v₀)")
					println("HELD Step 6 - Minimum Gibbs Energy = $(G₀)")
    			end
    			return  [1.0], [z₀], [v₀], G₀
    		end
				
			error = UBDⱽ - LBDⱽ 
			
			if verbose == true
        		println("HELD Step 3 - Update LBDⱽ from IPₓᵥ: LBDⱽ = $(LBDⱽ)")
        		println("HELD Step 3 - UBDⱽ - LBDⱽ = $(error) and tol = $(HELD_tol)")
    		end
			
			bphase = Vector{Float64}(undef,length(xmins_unique[1]))
			beta = Vector{Float64}(undef,length(xmins_unique))
			aphase = Matrix{Float64}(undef, length(xmins_unique[1]), length(xmins_unique))
			betaerror = 1.0
			if length(xmins_unique) > 1 
				for ib = 1:length(xmins_unique)
					sumx = 0.0
					for ia = 1:length(xmins_unique[1])-1
						aphase[ia,ib] = xmins_unique[ib][ia]
						sumx += aphase[ia,ib]
					end
					aphase[length(xmins_unique[1]),ib] = 1.0 - sumx
				end
				for ia = 1:length(xmins_unique[1])
					bphase[ia] = z₀[ia]
				end
				beta = aphase\bphase
				sumbeta = 0.0
				for ib = 1:length(xmins_unique)
					sumbeta += beta[ib]
				end
				betaerror = abs(1.0 - sumbeta)
				if verbose == true
        			println("HELD Step 3 - Test phase mole balance: error = $(betaerror) and tol =  $(sqrt(HELD_tol))")
        			println("HELD Step 3 - Phases found: np = $(length(beta))")
    			end
			end
			
			if verbose == true
        		println("HELD Step 3 - Test overall convergence:")
    		end
			
			np = length(beta)
			if error < HELD_tol && betaerror < sqrt(HELD_tol)
				if verbose == true
				    println("HELD Step 4 - Error within tolerences on UBDⱽ - LBDⱽ and phase mole balance")
        			println("HELD Step 4 - solution accepted")
    			end  			
	   			# normalise the solution before we do the Gibbs minimisation step.
	   			# its essential that xHELD moles balances and is a feasible solution
	   			phasemoles = Vector{Vector{Float64}}(undef,0)
	   			for ic = 1:nc-1			
					summoles = 0.0
					for ip = 1:np
						summoles += xmins_unique[ip][ic] * beta[ip]
					end
					mole = Vector{Float64}(undef,0)
					for ip = 1:np
						push!(mole, xmins_unique[ip][ic] * beta[ip] * z₀[ic] / summoles)
					end
					push!(phasemoles, mole)
				end
				
				summoles = 0.0
				for ip = 1:np
					x_nc = 1.0 - sum(xmins_unique[ip][1:nc-1])
					summoles += x_nc * beta[ip]
				end
				mole = Vector{Float64}(undef,0)
				for ip = 1:np
				    x_nc = 1.0 - sum(xmins_unique[ip][1:nc-1])
					push!(mole, x_nc * beta[ip] * z₀[nc] / summoles)
				end
				push!(phasemoles, mole)

				for ip = 1:np
					beta[ip] = 0.0
					for ic = 1:nc
						beta[ip] += phasemoles[ic][ip]
					end
				end

				for ip = 1:np
					for ic = 1:nc-1
						xmins_unique[ip][ic] = phasemoles[ic][ip] / beta[ip];
					end
				end
				# end of normalisation
				
	   			for ip = 1:np-1
    				push!(xHELD,beta[ip])
    				for ic = 1:nc
    					push!(xHELD,xmins_unique[ip][ic])
    				end
    			end
    			push!(xHELD,xmins_unique[np][nc])
    			HELD_complete = true
				break
			end
    		
    		if verbose == true
				println("HELD Step 3 - Add new (x,V)s: ℒs to the ℳ set and all current minimums to the ℳguess set")
    		end

			if error > 0.001
				limit_λs_by_bounds = true
			else
				limit_λs_by_bounds = false
				if verbose == true
					println("HELD Step 3 - λ bounds removed from OPₓᵥ")
				end
			end

			use_only_lowest_min = false
			if error > 0.001 && !use_only_lowest_min
				# add latest minimums to ℳguess
				for i = 1:length(fmins_unique)
					xminsGi_unique = append!(deepcopy(xmins_unique[i]),Gi(xmins_unique[i]))
					push!(ℳ,xminsGi_unique)
				end
			else
				# add lowest minimum to ℳ set
				for i = 1:length(ℒ)
					ℒGi = append!(deepcopy(ℒ[i]),Gi(ℒ[i]))
					push!(ℳ,ℒGi)
				end
			end

			# add minimums to ℳguess set
			# remove previous minimums from ℳguess only if number found is less than or equal to previous
			# this keeps the maximum number found in the set
			if length(fmins_unique) <= n_unique_previous
    			deleteat!(ℳguess, (length(ℳguess) - (n_unique_previous-1)):length(ℳguess))
    		end
    		# add latest minimums to ℳguess
    		for i = 1:length(fmins_unique)
				xminsGi_unique = append!(deepcopy(xmins_unique[i]),Gi(xmins_unique[i]))
				push!(ℳguess,xminsGi_unique)
			end

			n_unique_previous = length(fmins_unique)
			
			if verbose == true
        		println("HELD Step 3 - Overall convergence not satisfied return to step 2")
    		end
			   	
    	end
    	
    	if !HELD_complete
    	# we made it here without a HELD solution just return single phase and warn 	
    		if verbose == true
				println("HELD Step 5 - No solutions found")
				println("HELD Step 5 - Warning: fluid is being flagged as stable: try increasing max HELD iterations if this is a large component set problem")
				println("HELD Step 5 - Phase found 1")
				println("HELD Step 5 - Phase moles:")
				println("HELD Step 5 - Phase beta(1) = 1")
				println("HELD Step 5 - Phase mole fraction:")
				println("HELD Step 5 - Phase x(1) = $(z₀)")
				println("HELD Step 5 - Phase volumes:")
				println("HELD Step 5 - Phase volume(1) = $(v₀)")
				println("HELD Step 5 - Minimum Gibbs Energy = $(G₀)")
    		end
    		return  [1.0], [z₀], [v₀], G₀
    	else
    		# we made it here with a HELD solution use Gibbs Energy Minimisation to polish the solution this normally takes 2 to 3 iterations as we are close to the solution.
			if verbose == true
				println("HELD Step 5 - Start Gibbs Energy Minimisation:")
			end
			
			# xsol contains all the phases xgibbs works on np-1 phases and completes the missing phase via the mole balance.
			# to be sure we must make xgibbs feasible, no negatives and no greater than 1 values.
			
			Gibbs(x) = Gibbs_func(model,p,T,z₀,vref,np,x)
			Gibbs_g(x) = Solvers.gradient(Gibbs, x)
			Gibbs_h(x) = Solvers.hessian(Gibbs, x)

			lb = Vector{Float64}(undef,0)
		   	for ip = 1:np-1
				push!(lb,2.2e-16)
				for ic = 1:nc-1
					push!(lb,2.2e-16)
				end
				push!(lb,2.2e-16)
			end
			push!(lb,2.2e-16)
			
			ub = Vector{Float64}(undef,0)
		   	for ip = 1:np-1
				push!(ub,1.0 - 2.2e-16)
				for ic = 1:nc-1
					push!(ub,1.0 - 2.2e-16)
				end
				push!(ub,100.0)
			end
			push!(ub,100.0)
			
			projGibbs(x) = Projection(x,lb,ub)
			cnstGibbs(x,s) = Constraints(x,lb,ub,s)
			
			xsol,Gsol,iter,error,check = Solvers.trustregion_Dennis_Schnabel(Gibbs, Gibbs_g, Gibbs_h, projGibbs, cnstGibbs, xHELD, lb, ub, max_trust_region_iters, tol,false)
			
			if verbose == true
				println("HELD Step 5 - Gibbs Energy Minimisation: iterations taken = $(iter)")
				println("HELD Step 5 - Gibbs Energy Minimisation: error = $(error) - tol = $(tol)")
				if check
					println("HELD Step 5 - Gibbs Energy Minimisation: did not converge to required tolerance")
					if error < HELD_tol && iter < max_trust_region_iters
						println("HELD Step 5 - Gibbs Energy Minimisation: solution accepted as this is less than HELD tolerence $(HELD_tol)")
					else
						println("HELD Step 5 - Gibbs Energy Minimisation: not converged")
					end
				else
					println("HELD Step 5 - Gibbs Energy Minimisation: solution found")
				end
			end
			
			# unpack xsol using mole balance
			beta = Vector{Float64}(undef,0)
			xp = Vector{Vector{Float64}}(undef,0)
			vp = Vector{Float64}(undef,0)
			ix = 0
			for ip=1:np-1
				ix += 1
				push!(beta,xsol[ix])
				xc = Vector{Float64}(undef,0)
				for ic = 1:nc-1
					ix += 1
					push!(xc,xsol[ix])
				end
				push!(xc,1.0 - sum(xc[1:nc-1]))
				push!(xp,xc)
				ix += 1
				push!(vp,vref/xsol[ix])
			end
			ix += 1
			# volume is now returned
			push!(vp,vref/xsol[ix])
			push!(beta,1.0 - sum(beta[1:np-1]))
			xc = Vector{Float64}(undef,0)
			for ic = 1:nc
				sum_moles = 0.0
				for ip = 1:np-1
					sum_moles += beta[ip]*xp[ip][ic]
				end
				push!(xc,(z₀[ic] - sum_moles)/beta[np])
			end
			push!(xp,xc)
			
			if verbose == true
				println("HELD Step 5 - Gibbs Energy Minimisation: Normalise final solution so it mole balances")
			end
			# normalise the solution before we finish.
		   	phasemoles = Vector{Vector{Float64}}(undef,0)
		   	for ic = 1:nc			
				summoles = 0.0
				for ip = 1:np
					summoles += xp[ip][ic] * beta[ip]
				end
				mole = Vector{Float64}(undef,0)
				for ip = 1:np
					push!(mole, xp[ip][ic] * beta[ip] * z₀[ic] / summoles)
				end
				push!(phasemoles, mole)
			end

			for ip = 1:np
				beta[ip] = 0.0
				for ic = 1:nc
					beta[ip] += phasemoles[ic][ip]
				end
			end

			for ip = 1:np
				for ic = 1:nc-1
					xp[ip][ic] = phasemoles[ic][ip] / beta[ip];
				end
				ρp = HELD_density(model,p,T,xp[ip],vref)
				vp[ip] = vref/ρp
			end
			# end of 
			
			ivpsort = sortperm(vp, rev = true)

			betas = Vector{Float64}(undef,0)
			for ip in eachindex(beta)
				push!(betas,beta[ivpsort[ip]])
			end

			xps = Vector{Vector{Float64}}(undef,0)
			for ip in eachindex(xp)
				# undo sort
				xs = Vector{Float64}(undef,length(xp[ivpsort[ip]]))
				for i in eachindex(xp[ivpsort[ip]])
					xs[i] = xp[ivpsort[ip]][i]
				end 
				push!(xps,xs)
			end

			vps = Vector{Float64}(undef,0)
			for ip in eachindex(vp)
				push!(vps,vp[ivpsort[ip]])
			end
			
			if verbose == true
				println("HELD Step 5 - Complete")
				println("HELD Step 6 - Phases found $(length(betas))")
				println("HELD Step 6 - Phase moles:")
				for ip = 1:length(betas)
					println("HELD Step 6 - Phase beta[$(ip)] = $(betas[ip])")
				end
				println("HELD Step 6 - Phase mole fraction:")
				for ip = 1:length(betas)
					println("HELD Step 6 - Phase x[$(ip)] = $(xps[ip])")
				end
				println("HELD Step 6 - Phase volumes:")
				for ip = 1:length(betas)
					println("HELD Step 6 - Phase volume[$(ip)] = $(vps[ip])")
				end
				println("HELD Step 6 - Minimum Gibbs Energy = $(Gsol)")
			end
			
			return betas,xps,vps,Gsol
			
		end # Gibss Minimisation
		
    end # Search for phases
    
end

function HELD_clean_local_solutions(G₀, x₀, fmins, xmins, tol, verbose)
    iremove = fill(false,length(fmins))
    iminfound = fill(false,length(fmins))
    for ir = 1:length(fmins)    	
    	if !iremove[ir] && !iminfound[ir]
    	    imin = ir
    		fmin = fmins[ir]
    		for imins = 1:length(fmins)
    			if !iremove[imins] && !iminfound[imins]
        			if fmins[imins] < fmin
        				imin = imins
        			    fmin = fmins[imins]
    				end
    			end
    		end
			iminfound[imin] = true
    		for imins = 1:length(fmins)
    			if !iremove[imins] && !iminfound[imins]
    				distances = xmins[imins] - xmins[imin]
    				distance = norm(distances, Inf)
    				if distance < sqrt(tol)
    					iremove[imins] = true
    				end
    			end
    		end
    	end
    end
    # remove trival solutions
    for imins = 1:length(fmins)
    	if iminfound[imins]
    		distances = xmins[imins] - x₀
    		distance = norm(distances, Inf)
    		if distance < sqrt(tol)
    			iminfound[imins] = false
    		end
    	end
    end
    fmins_unique = Vector{Float64}(undef,0)
    xmins_unique = Vector{Vector{Float64}}(undef,0)
    stable = true
    for ir = 1:length(fmins)
    	if iminfound[ir]
    		if fmins[ir] < G₀
    			stable = false
    		end 
    		push!(fmins_unique,fmins[ir])
    		push!(xmins_unique,xmins[ir])
    	end
    end
    return fmins_unique, xmins_unique, stable
end

export HELDTPFlash
