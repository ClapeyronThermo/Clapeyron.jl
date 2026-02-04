

function ADScalarObjective(f,x0::T) where T <: AbstractArray
    Hresult = DiffResults.HessianResult(x0)
    Hconfig = ForwardDiff.HessianConfig(f,Hresult,x0)
    function fg(_df,x)
        Gconfig = _GradientConfig(Hconfig)
        result = ForwardDiff.gradient!(Hresult,f,x,Gconfig)
        if __is_implace(x)
            _df .= DiffResults.gradient(result)
            df = _df
        else
            df = DiffResults.gradient(result)
        end
        fx = DiffResults.value(result)
        return fx,df
    end

      function g(df,x)
        _,df = fg(df,x)
        return df
    end

    function fgh(df,d2f,x)
        result = ForwardDiff.hessian!(Hresult,f,x,Hconfig)
        fx = DiffResults.value(result)
        if __is_implace(x)
            d2f .= DiffResults.hessian(result)
            df .= DiffResults.gradient(result)
            return fx,df,d2f
        else
            _d2f = DiffResults.hessian(result)
            _df = DiffResults.gradient(result)
            return fx,_df,_d2f
        end 
    end

    function h(d2f,x)
        result = ForwardDiff.hessian!(Hresult,f,x)
        if __is_implace(x)
            d2f .= DiffResults.hessian(result)
            return d2f
        else
            return DiffResults.hessian(result)
        end
    end

    return ScalarObjective(f=f,
    g=g,
    fg=fg,
    fgh=fgh,
    h=h)
end

function ADScalarObjective(f::F,x0::T) where {F,T <: SVector}
    function g(_a1,x)
        return ForwardDiff.gradient(f,x)
    end
    
    function fg(_a1,x)
        gres = DiffResults.GradientResult(x)
        result = ForwardDiff.gradient!(gres,f,x)
        df = DiffResults.gradient(result)
        fx = DiffResults.value(result)
        return fx,df
    end

    function fgh(_a1,_a2,x)
        Hresult = DiffResults.HessianResult(x)
        result = static_fgh(Hresult,f,x)
        fx = DiffResults.value(result)
        _d2f = DiffResults.hessian(result)
        _df = DiffResults.gradient(result)
        return fx,_df,_d2f
    end

    function h(_a1,x)
        return ForwardDiff.hessian(f,x)
    end

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

    function fg(_,x)
        ∂fx,x = f∂f(f,x)
        return ∂fx,x
    end

    function fgh(_,_,x)
        fx,∂fx,∂2fx = f∂f∂2f(f,x)
        return fx,∂fx,∂2fx
    end

    function h(_,x)
        ∂2fx = derivative(g,x)
        return ∂2fx
    end

    return ScalarObjective(f = f,
    g = g,
    fg = fg,
    fgh = fgh,
    h = h)
end

#uses brent, the same default that Optim.jl uses
function optimize(f,x0::NTuple{2,T},method=BrentMin(),options=OptimizationOptions()) where T<:Real
    scalarobj = ADScalarObjective(f,first(x0))
    optprob = OptimizationProblem(obj = scalarobj,bounds = x0, inplace=false)
    res = NLSolvers.solve(optprob,method,options)
    return res
end

#general one, with support for ActiveBox
function optimize(f,x0,method=LineSearch(Newton2(x0),NLSolvers.Static(1.0)),options=OptimizationOptions();bounds = nothing)
    inplace = __is_implace(x0)
    scalarobj = ADScalarObjective(f,x0)
    optprob = OptimizationProblem(obj = scalarobj,inplace = inplace,bounds = bounds)
    return NLSolvers.solve(optprob,x0,method,options)
end

function optimize(optprob::OptimizationProblem,x0,method=LineSearch(Newton2(x0),NLSolvers.Static(1.0)),options=OptimizationOptions();bounds = nothing)
    return NLSolvers.solve(optprob,x0,method,options)
end

#build scalar objective -> Optimization Problem
function optimize(scalarobj::ScalarObjective,x0,method=LineSearch(Newton2(x0),NLSolvers.Static(1.0)),options=OptimizationOptions();bounds = nothing)
    inplace = __is_implace(x0)
    optprob = OptimizationProblem(obj = scalarobj,inplace = inplace,bounds = bounds)
    return NLSolvers.solve(optprob,x0,method,options)
end

function optimize(f,x0,method::NLSolvers.NelderMead,options=OptimizationOptions();bounds = nothing)
    scalarobj = ScalarObjective(f = f)
    optprob = OptimizationProblem(obj = scalarobj,inplace = false,bounds = bounds)
    return NLSolvers.solve(optprob,x0,method,options)
end


#= only_fg!: Optim.jl legacy form:
function fg!(F,G,x)
  # do common computations here
  # ...
  if G != nothing
    # code to compute gradient here
    # writing the result to the vector G
  end
  if F != nothing
    # value = ... code to compute objective function
    return value
  end
end
=#

function only_fg!(fg!::T) where T
    function f(x)
        return fg!(true,nothing,x)
    end

    function g(df,x)
        fg!(nothing,df,x)
        return df
    end
    function fg(df,x)
        fx = fg!(true,df,x)
        return fx,df
    end

    return ScalarObjective(f=f,
    g=g,
    fg=fg,
    fgh=nothing,
    h=nothing)
end

#= only_fgh!: Optim.jl legacy form:
function fgh!(F,G,H,x)
  G == nothing || # compute gradient and store in G
  H == nothing || # compute Hessian and store in H
  F == nothing || return f(x)
  nothing
end
=#

function only_fgh!(fgh!::T) where T
    function f(x)
        return fgh!(true,nothing,nothing,x)
    end

    function g(df,x)
        fgh!(nothing,df,nothing,x)
        return df
    end

    function fg(df,x)
        fx = fgh!(true,df,nothing,x)
        return fx,df
    end

    function fgh(df,d2f,x)
        fx = fgh!(true,df,d2f,x)
        return fx,df,d2f
    end

    function h(d2f,x)
        fgh!(nothing,nothing,d2f,x)
        return d2f
    end

    return ScalarObjective(f=f,
    g=g,
    fg=fg,
    fgh=fgh,
    h=h)
end

struct BoundOptim1Var end

function optimize(f,x0::NTuple{2,T},method::BoundOptim1Var,options=OptimizationOptions()) where T<:Real
    return _1var_optimize_quad(f,x0,options)
end

quad_interp(x,f) = quad_interp(x[1],x[2],x[3],f[1],f[2],f[3])

function quad_interp(xa,xb,xc,fa,fb,fc)
        #f1 = ax12 + bx1 + c
    #f2 = ax22 + bx2 + c
    #f3 = ax32 + bx3 + c
    A = @SMatrix [xa*xa xa oneunit(xa); xb*xb xb oneunit(xb); xc*xc xc oneunit(xc)]
    B = SVector((fa,fb,fc))
    z = A\B
    a,b,c = z
    return a,b,c
end

function mean2(xhi,xlo)
    dx = xhi - xlo
    return xlo + 0.6dx
end

function _1var_optimize_quad(f,x0,options)
    xlo, xhi = minmax(x0[1],x0[2])
    time0 = time()
    # Initialize three points
    x1 = xlo
    x3 = xhi
    x2 = (x1 + x3) / 2
    xlo0 = xlo
    xhi0 = xhi
    f1 = f(x1)
    f2 = f(x2)
    f3 = f(x3)
    f30 = f3
    fx = NaN*f1
    xx = NaN*f1
    iter_x = 0
    not_improved = 0
    for iter in 1:options.maxiter
        iter_x += 1
        iter_x == options.maxiter && break
        not_improved == 2 && break
        a,b,c = quad_interp(x1,x2,x3,f1,f2,f3)
        xlo,xhi = extrema((x1, x2, x3))
        df = f30 - fx
        if a <= 0 || abs(xlo - xhi) < options.x_abstol
            # Parabola opens downward or is linear, use best current point
            fvals = SVector((f1, f2, f3))
            xvals = (x1, x2, x3)
            f_best,idx = findmin(fvals)
            minidx = argmin(fvals)
            x_best = xvals[idx]
            xx = x_best
            fx = f_best
            break
        end
        x_interp = -b/2a
        
        # Ensure new point is within current bounds
        x_new = xlo0 <= x_interp <= xhi0 ? x_interp : mean2(xhi,xlo)
        x_new = clamp(x_interp,xlo0,xhi0)
        f_new = f(x_new)

        # Create array of all four points

        all_f = SVector((f1,f2,f3,f_new))
        all_x = (x1,x2,x3,x_new)
        idxs = sortperm(all_f)
        i1,i2,i3 = idxs[1],idxs[2],idxs[3]
        # Keep the three best points (reject worst)
        x1, f1 = all_x[i1],all_f[i1]
        x2, f2 = all_x[i2],all_f[i2]
        x3, f3 = all_x[i3],all_f[i3]
        if xx == x1
            not_improved += 1
        end
        xx = x1
        fx = f1
    end
    
    return NLSolvers.ConvergenceInfo(
        BoundOptim1Var(),
        (; solution = xx, f0 = f3, minimum = fx, time = time() - time0, iter = iter_x),
        options,
    )
end
