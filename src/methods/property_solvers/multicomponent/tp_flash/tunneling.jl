function rand_tunneling!(method::Nothing,x)
    x0 = rand(length(x))
    ∑x = sum(x0)
    x0 .= x0./∑x
    x0[1] = rand() #to break the correlation between all the values and the volume
    return x0
end

function eos_tunneling(model,p,T,n,verbose = false,p̃ = 1e-3)
    nc = length(model)
    # Stage I: Stability test and initialisation
    if verbose
        println("==========================================")
        println("Stage I: Stability test and initialisation")
        println("==========================================")
        println("----------------------------")
        println("Step 1: Stability test at n₀")
        println("----------------------------")
    end
    nn = n/sum(n)
    vₛ = volume(model,p,T,nn)
    μ₀ = VT_chemical_potential(model,vₛ,T,nn)
    xₛ = deepcopy(nn)
    xₛ[1] = lb_volume(model,nn)/vₛ

    d(x) = tunneling_d(model,T,x,μ₀)
    if verbose
        println("Initial point found. Beginning tunneling")
    end
    d0 = d(xₛ)
    f = TunnelingF(d,d0,deepcopy(xₛ),p̃)
    x0 = zeros(typeof(d0),length(μ₀))
    for _ in 1:10nc
        x0 = rand_tunneling!(nothing,x0) #TODO: add a deterministic sampling method
        r = Solvers.optimize(f,x0)
        xₙ = Solvers.x_sol(r)
        dₙ = d(xₙ)
        if dₙ<0
            verbose && println("Negative tangent found. Moving on to step 2.")

            return (xˢ,gibbs_free_energy(model,p,T,xˢ))
        elseif dₙ < f.fstar #lower known minimum found
            f.xstar .= xₙ
            f.fstar = dₙ
        end
    end
    if verbose
        println("No negative tangent found. Initial point is stable.")
        println("Terminating HELD")
    end
    return (xˢ,gibbs_free_energy(model,p,T,xˢ))
end

mutable struct TunnelingF{F,T,X}
    f::F
    fstar::T
    xstar::X
    ptilde::Float64
end

function (state::TunnelingF)(x)
    f = state.f(x)
    fstar = state.fstar
    xstar = state.xstar
    p̃ = state.ptilde
    normx = √(sum((x[i]-xstar[i] for i in eachindex(x))^2))
    return (f(x)-fstar)*exp(p̃/normx)
end

function tunneling_d(model,T,y,μ₀)
    x = FractionVector(@view(y[2:end]))
    lb_v = lb_volume(model,x)
    V = (lb_v/y[1])
    ∑x = sum(x)
    ∑xᵢμ₀ᵢ = dot(x,μ₀)
    return (eos(model,V,T,x)+p*V)/∑x - ∑xᵢμ₀ᵢ
end


