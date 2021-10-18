
#=
struct TunnelingModifier{T,F,V}
    f::F
    x::V
    fbest::Ref{T}
end


function tunneling(f,lb,ub,x0)
    N = length(ub)
    # Relevant configuration
    tolf=1e-8
    options = OptimizationOptions(g_abstol=5e-8)
    prob = ADScalarObjective(f,x0)
    # Minimisation phase
    prob = OptimizationProblem(obj=scalarobj, 
    bounds=(lb,ub); inplace=true)
    opt_min.xtol_rel = 1e-8
    res =  NLSolvers.solve(prob, x0, ModifiedActiveBox(Newton();epsilon=tolbounds), options)
    

    best_f = res.info.minimum
    opt_x  = []
    best_x = res.info.minimizer
    append!(opt_x,[min_x])

    # Tunneling phase
    opt_tun = NLopt.Opt(:LD_MMA, length(ub))
    opt_tun.lower_bounds = lb
    opt_tun.upper_bounds = ub
    opt_tun.xtol_rel = 1e-8
    opt_tun.stopval = -1e-6

    for i in 1:10*N
        T0 = x -> (f(x)-f_best)*prod(exp(1e-1/sqrt(sum((x[i]-x_opt[j][i])^2 for i in 1:N))) for j in 1:k)
        T  = (x,g) -> NLopt_obj(T0,x,g)
        opt_tun.min_objective =  T

        r  = 2.0.*rand(Float64,(N)).-1.0
        ϵ1 = 2*(tolf)^(1/5)*(1+norm(best_x,2))
        x0 = r/norm(r,2).*ϵ1+best_x
        x0 = ub.*(x0.>=ub)+lb.*(x0.<=lb)+x0.*(ub.>x0.>lb)

        # Tunneling
        (new_f,new_x,status) = NLopt.optimize(opt_tun, x0)
        if minimum(res.info.∇fz) <=
            println(i)
            break
        end
        # Minimisation
        (min_f,min_x,status) = NLopt.optimize(opt_min, new_x)
        if min_f<best_f
            best_f = min_f
            best_x = min_x
            opt_x  = []
            append!(opt_x,[min_x])
        else
            append!(opt_x,[min_x])
        end
    end
    return (best_f,best_x)
end

function NLopt_obj(f,x,g)
        if length(g) > 0
            df = DiffResults.GradientResult(x)
            df = ForwardDiff.gradient!(df,f,x)
            g .= DiffResults.gradient(df)
            return DiffResults.value(df)
        else
            return f(x)
        end
end

=#