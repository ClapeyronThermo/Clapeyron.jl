
function ADScalarObjective(f,x0::AbstractArray)
    Hres = DiffResults.HessianResult(x0)
    function _g(df,x,Hresult)
        ForwardDiff.gradient!(Hresult,f,x)
        df .= DiffResults.gradient(Hresult)
        df
    end
    function _fg(df,x,Hresult)
        ForwardDiff.gradient!(Hresult,f,x)
        df .= DiffResults.gradient(Hresult)
        fx = DiffResults.value(Hresult)
        return fx,df
    end

    function _fgh(df,d2f,x,Hresult)
        ForwardDiff.hessian!(Hresult,f,x)
        d2f .= DiffResults.hessian(Hresult)
        df .= DiffResults.gradient(Hresult)
        fx = DiffResults.value(Hresult)
        return fx,df,d2f
    end

    function h(d2f,x)
        ForwardDiff.hessian!(d2f,f,x)
        d2f
    end
    g(df,x) = _g(df,x,Hres)
    fg(df,x) = _fg(df,x,Hres)
    fgh(df,d2f,x) =  _fgh(df,d2f,x,Hres)
    return ScalarObjective(f=f,
    g=g,
    fg=fg,
    fgh=fgh,
    h=h)
end

function ADScalarObjective(f,x0::Number)
    function g(x)
        return derivative(f,x)
    end

    function fg(∂fx,x)
        ∂fx,x = f∂f(f,x)
        return ∂fx,x
    end

    function fgh(∂fx,∂2fx,x)
        fx,∂fx,∂2fx = f∂f∂2f(f,x)
        return fx,∂fx,∂2fx
    end

    function h(∂2fx,x)
        ∂2fx = derivative(g,x)
        return ∂2fx
    end

    return ScalarObjective(f=f,
    g=g,
    fg=fg,
    fgh=fgh,
    h=h)
end
"""
    function optimize(f,x0,method=LineSearch(Newton()), options=OptimizationOptions())
"""
function optimize(f,x0::Number,method=LineSearch(Newton()),options=OptimizationOptions())
    scalarobj = ADScalarObjective(f,x0)   
    optprob = OptimizationProblem(scalarobj; inplace=false) 
    return NLSolvers.solve(optprob, x0, method,options)
end

#this should allow to use another schemes to update
struct ModifiedActiveBox{M,T}
    scheme::M
    ϵ::T
end

Base.summary(ab::ModifiedActiveBox) = "Active Box, using $(Base.summary(ab.scheme))"

NLSolvers.modelscheme(ab::ModifiedActiveBox) = ab.scheme
ModifiedActiveBox(scheme=Newton(); epsilon=1e-12) =  ModifiedActiveBox(scheme,epsilon)

function box_optimize(f,x0,lb,ub,options=OptimizationOptions();tolbounds = 1e-12)
    scalarobj = ADScalarObjective(f,x0)
    prob = OptimizationProblem(obj=scalarobj, 
    bounds=(lb,ub); inplace=true)
    NLSolvers.solve(prob, x0, ModifiedActiveBox(Newton(linsolve = cholesky_linsolve);epsilon=tolbounds), options)
end

#see https://github.com/JuliaNLSolvers/NLSolvers.jl/issues/21
function NLSolvers.solve(prob::OptimizationProblem, x0, scheme::ModifiedActiveBox, options::OptimizationOptions)
    t0 = time()

    x0, B0 = x0, (false*x0*x0' +I)
    lower, upper = NLSolvers.bounds(prob)
    ϵbounds = mapreduce(b->(b[2] - b[1])/2, min, zip(lower, upper)) # [1, pp. 100: 5.41]

    !any(clamp.(x0, lower, upper) .!= x0) || error("Initial guess not in the feasible region")

    linesearch = NLSolvers.ArmijoBertsekas()
    mstyle = NLSolvers.OutOfPlace()
    direction_scheme = NLSolvers.modelscheme(scheme)
    objvars = NLSolvers.prepare_variables(prob, scheme, x0, copy(x0), B0)
    f0, ∇f0 = objvars.fz, norm(objvars.∇fz, Inf) # use user norm
    fz, ∇fz = objvars.fz, objvars.∇fz # use user norm
    fx, ∇fx = fz, copy(∇fz)
    B = B0
    x, z = copy(x0), copy(x0)
    Tf = typeof(fz)
    is_first=false
    Ix = Diagonal(z.*0 .+ 1)

    unconstrained_d = copy(∇fz)
    d = copy(∇fz)
    dmin = copy(∇fz)
    y = nothing
    s = copy(∇fz)
    xnorm = copy(∇fz)
    activeset = similar(∇fz,Bool)
    Ĥ = copy(B)
    Ĥchol = cholesky(Ĥ)
    for iter = 1:options.maxiter
        x .= z
        fx = fz
        ∇fx .= ∇fz
        
        ϵ = min(norm(clamp.(x.-∇fx, lower, upper).-x), ϵbounds) # Kelley 5.41 and just after (83) in [1]
        activeset .= NLSolvers.is_ϵ_active.(x, lower, upper, ∇fx, ϵ)
        Ĥ .= NLSolvers.diagrestrict.(B, activeset, activeset', Ix)
        #d = clamp.(x.-HhatChol\∇fx, lower, upper).-x
        # Update current gradient and calculate the search direction
        NLSolvers.find_direction!(unconstrained_d,Ĥ,objvars.Pg,∇fx,direction_scheme)
        d .= clamp.(x .+ unconstrained_d, lower, upper).-x #clamp unconstrained direction
        dmin .= d./unconstrained_d
        d_scale = minimum(dmin) 
        #scale the direction, the original direction is the same.
        #d .=  @. unconstrained_d*d_scale + (one(d_scale)-d_scale)*d
        !iszero(d_scale) && (d .= unconstrained_d .* d_scale)
        #if the scale is zero, we are actually in a boundary, and we can only grow parallel along the box.
        φ = NLSolvers._lineobjective(mstyle, prob, ∇fz, z, x, d, fz, dot(∇fz, d))
        
        # Perform line search along d
        # Also returns final step vector and update the state
        α, f_α, ls_success = NLSolvers.find_steplength(mstyle, linesearch, φ, Tf(1), ∇fz, activeset, lower, upper, x, d, ∇fx, activeset)
        # # Calculate final step vector and update the state
        if !ls_success
            return NLSolvers.ConvergenceInfo(scheme, (prob=prob, B=B, ρs=norm(x.-z), ρx=norm(x), minimizer=z, fx=fx, minimum=fz, ∇fz=∇fz, f0=f0, ∇f0=∇f0, iter=iter, time=time()-t0), options)
        end

        s .= @. α * d
        z .= @. x + s
        s .= clamp.(z, lower, upper) .- x
        z .= x .+ s
        # Update approximation. it only works with newton at the moment
        fz, ∇fz, B, s, y = NLSolvers.update_obj!(prob.objective,s,y, ∇fx, z, ∇fz, B, direction_scheme)
        #fz, ∇fz, B, s, y = NLSolvers.update_obj(prob.objective, s, ∇fx, z, ∇fz, B, direction_scheme,is_first)
        #@show fz, ∇fz, B, s, y
        #@show Ĥ
        xnorm .= x.-clamp.(x.-∇fz, lower, upper) 
        normx = norm(xnorm, Inf)
        if normx < options.g_abstol
            return NLSolvers.ConvergenceInfo(scheme, (prob=prob, B=B, ρs=norm(x.-z), ρx=norm(x), minimizer=z, fx=fx, minimum=fz, ∇fz=∇fz, f0=f0, ∇f0=∇f0, iter=iter, time=time()-t0), options)
        end
    end
  iter =  options.maxiter
  return NLSolvers.ConvergenceInfo(scheme, (prob=prob, B=B, ρs=norm(x.-z), ρx=norm(x), minimizer=z, fx=fx, minimum=fz, ∇fz=∇fz, f0=f0, ∇f0=∇f0, iter=iter, time=time()-t0), options)
end


#=
julia> Optim.minimizer(results)
2-element Vector{Float64}:
 1.2500000000000002
 1.5626548683420372

julia> Optim.minimum(results)
0.06250239842033663
=#