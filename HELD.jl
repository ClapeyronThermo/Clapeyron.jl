using Optim


function Rosenbrock(x::Vector)
    N=10
    f(x) = sum(100*(x[i+1]-x[i]^2)^2+(1-x[i])^2 for i in 1:N-1)
    return f(x)
end

function Rastrigin(x::Vector)
    N=2
    f(x) = 10*N.+sum(x[i]^2-10*cos(2*pi*x[i]) for i in 1:N)
    return f(x)
end

function Schubert(x::Vector)
    f(x) = sum(i*cos((i+1)*x[1]+i) for i in 1:5)*sum(i*cos((i+1)*x[2]+i) for i in 1:5)
    return f(x)
end

function Tunneling(obj_f,lb,ub)
    N = length(ub)
    # Relevant configuration
    # opt = NLopt.Opt(:LD_MMA, length(ub))
    # opt.lower_bounds = lb
    # opt.upper_bounds = ub
    # opt.ftol_abs = 1e-8
    # opt.stopval = -1e-8

    # Minimisation phase
    # opt.min_objective =  obj_f
    inner_optimizer = GradientDescent()
    x0 = (ub.-lb).*rand(Float64,(N)).+lb
    res =optimize(obj_f, lb, ub, x0, Fminbox(inner_optimizer); autodiff = :forward)
    f_best = Optim.minimum(res)
    x_best = Optim.minimizer(res)
    x_opt = []
    append!(x_opt,[x_best])
    println("Minimisation complete")
    # Tunneling phase
    for k in 1:5

        T(x) = (obj_f(x)-f_best)*prod(exp(1e-1/sqrt(sum((x[i]-x_opt[j][i])^2 for i in 1:N))) for j in 1:k)

        x0 = (2.0.*rand(Float64,(N)).-1.0).*1e-2+x_best
        res =optimize(T, lb, ub, x0, Fminbox(inner_optimizer); autodiff = :forward)
        f_new = Optim.minimum(res)
        x_new = Optim.minimizer(res)
        append!(x_opt,[x_new])
        if f_new<0
            f_best = obj_f(x_new)
            x_best = x_new
        end
        println("Tunneling stage complete")
    end
    return (x_best,f_best)
end
N=2
# (min_f,min_x,status) = Tunneling(Rosenbrock,-2.048*ones(N),2.048*ones(N))
(x,f) = Tunneling(Schubert,-10*ones(N),10*ones(N))

println(x,f)
