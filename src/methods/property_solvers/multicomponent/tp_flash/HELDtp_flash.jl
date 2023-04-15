using JuMP, HiGHS

"""
    HELDTPFlash(;numphases = 2;max_steps = 1e4*(numphases-1),population_size =20,time_limit = Inf,verbose = false, logspace = false)

Method to solve non-reactive multicomponent flash problem by finding global minimum of Gibbs Free Energy via Differential Evolution.

User must assume a number of phases, `numphases`. If true number of phases is smaller than numphases, model should predict either (a) identical composition in two or more phases, or (b) one phase with negligible total number of moles. If true number of phases is larger than numphases, a thermodynamically unstable solution will be predicted.

The optimizer will stop at `max_steps` evaluations or at `time_limit` seconds

"""
Base.@kwdef struct HELDTPFlash <: TPFlashMethod
    outer_steps::Int = 100
    inner_steps::Int = 100
    time_limit::Float64 = Inf
    eps_λ::Float64 = 0.5
    eps_b::Float64 = 1e-2
    eps_x::Float64 = 1e-3
    eps_η::Float64 = 1e-3
    eps_μ::Float64 = 1e-6
    eps_g::Float64 = 1e-6
    p_tun::Float64 = 1e-3
    verbose::Bool = false
    logspace::Bool = false
end

index_reduction(flash::HELDTPFlash,z) = flash

#include("tunneling.jl")

function _Lⱽ(model,p,T,x,λˢ,n) #x = [log10(v),w]
    V = exp10(x[1])
    w = Fractions.FractionVector(@view(x[2:end]))
    nc = length(model)
    (eos(model,V,T,w)+p*V)/R̄/T + sum(λˢ[j]*(n[j]-x[j+1]) for j ∈ 1:nc-1)
end

#Inner problem cache
struct HELDIPCache{𝕋}
    model::𝕋 #current eos model
    p::Float64 #pressure of the system, pa
    T::Float64 #temperature of the system, K
    z::Vector{Float64} #composition of the feed
    ℳ::Vector{Vector{Float64}} #list of phases in the form [log10(v);w]
    G::Vector{Float64} #gibbs energies of the phases
    λᴸ::Vector{Float64} #lower bound of lagrangians
    λᵁ::Vector{Float64} #upper bound of lagrangians
    μ::Vector{Vector{Float64}} #list of chemical potentials
    LV::Vector{Float64} #list of tpd for each phase
end

Base.length(cache::HELDIPCache) = length(cache.G)

#constructor
function HELDIPCache(model,p,T,z)
    nc = length(z)
    ℳ = Vector{Vector{Float64}}(undef,0)
    G = Vector{Float64}(undef,0)
    λᴸ = fill(Inf,nc - 1)
    λᵁ = fill(-Inf,nc - 1)
    μ = Vector{Vector{Float64}}(undef,0)
    LV = Vector{Float64}(undef,0)
    return HELDIPCache(model,p,T,z,ℳ,G,λᴸ,λᵁ,μ,LV)
end

function add_candidate!(cache::HELDIPCache,xx;precalc = false)
    p,T = cache.p,cache.T
    model = cache.model
    nc = length(model)

    if precalc
        v = exp10(first(xx))
        x = xx[2:end]
        m = xx
    else
        v = volume(model,p,T,xx)
        x = xx
        m = vcat(log10(v),x)
    end
    μ_new = VT_chemical_potential(model,v,T,x)/R̄/T
    
    G_new = dot(μ_new,x)
    if isnan(G_new)
        @show μ_new
        @show VT_gibbs_free_energy(model,v,T,x)/R̄/T
        @show x,v
    end
    μλ = copy(μ_new)
    μλ_nc = μλ[end]
    μλ .-= μλ_nc
    resize!(μλ,nc-1)
    nc = length(cache.model)
    #adding the calculated values:
    push!(cache.ℳ,m)
    push!(cache.G,G_new)
    update_bounds!(cache::HELDIPCache,μλ)
    push!(cache.μ,μ_new)
    push!(cache.LV,NaN)
    return cache
end

function already_in_cache(cache::HELDIPCache,xx,precalc = false)
    if !precalc
        v = volume(model,p,T,xx)
        x = xx
        m = vcat(log10(v),x)
    else
        m = xx
    end
    for xi in cache.ℳ
        if isapprox(xi,m,rtol = 1e-2)
            return true
        end
    end
    return false
end

function update_bounds!(cache::HELDIPCache,λ)
    nc = length(cache.model)
    for i in 1:nc-1
        λi = λ[i]
        if cache.λᴸ[i] > λi
            cache.λᴸ[i] = λi
        end
        if cache.λᵁ[i] < λi
            cache.λᵁ[i] = λi
        end
    end
    return cache
end

function add_candidate_with_λ!(cache::HELDIPCache,x,λ)
    @show λ
    cache = add_candidate!(cache,x)
    xx = cache.ℳ[end]
    model,p,T = cache.model,cache.p,cache.T
    n = cache.z
    cache.LV[end] = _Lⱽ(model,p,T,xx[1:end-1],λ,n)
    return cache,true
end

function tp_flash_impl(model::EoSModel, p, T, n, method::HELDTPFlash)
    nc = length(model)
    p̃ = method.p_tun
    TYPE = typeof(p+T+first(n))
    # Stage I: Stability test and initialisation
    if method.verbose == true
        println("==========================================")
        println("Stage I: Stability test and initialisation")
        println("==========================================")
        println("----------------------------")
        println("Step 1: Stability test at n₀")
        println("----------------------------")
    end
    v00 = volume(model,p,T,n)
    vₛ = log10(v00)
    μ₀ = VT_chemical_potential(model,v00,T,n)
    xₛ = prepend!(deepcopy(n),vₛ)

    d(x) = (eos(model,10 .^x[1],T,x[2:end])+p*10 .^x[1])/sum(x[2:end])-∑(x[2:end].*μ₀)
    if method.verbose == true
        println("Initial point found. Beginning tunneling")
    end
    i = 0
    while i<10nc
        f(x) = d(x)*exp(p̃/√(∑((xₛ[2:end]-x[2:end]).^2)))
        x0 = rand(length(n))
        x0 = x0./sum(x0)
        v0 = log10(volume(model,p,T,x0))
        x0 = prepend!(x0,v0)
        r = Solvers.optimize(f,x0)
        xₙ = Solvers.x_sol(r)
        dₙ = d(xₙ)
        if dₙ<0
            if method.verbose == true
                println("Negative tangent found. Moving on to step 2.")
            end
            break
        end
        i+=1
    end
    if i==10nc
        if method.verbose == true
            println("No negative tangent found. Initial point is stable.")
            println("Terminating HELD")
        end
        return (xˢ,gibbs_free_energy(model,p,T,xˢ))
    end

    if method.verbose == true
        println("--------------------------------------")
        println("Step 2: Initialisation of dual problem")
        println("--------------------------------------")
    end
    k = 0
    UBDⱽ = gibbs_free_energy(model,p,T,n)/R̄/T
    cache = initial_candidate_phases(model,p,T,n,μ₀)
    #λ₀ = (μ₀[1:nc-1].-μ₀[nc])/R̄/T
    #LV = vec(G.+ sum(λ₀'.*(n[1:nc-1]'.-ℳ[:,1:nc-1]),dims=2))
    if method.verbose == true
        println("Iteration counter set to k="*string(k))
        println("Upper bound set to UBDⱽ="*string(UBDⱽ))
        println("ℳ initialised")
        println("===================================================")
        println("Stage II: Identification of candidate stable phases")
        println("===================================================")
    end
    nps=1
    
    active_r = Ref{Vector{Bool}}()
    while k<=method.outer_steps
        UBDⱽ,cache,active = HELD_stage_II(cache,UBDⱽ,method,k)
        active_r[] = active
        nps = count(active)
        if nps>=2
            break
        end
        k+=1
    end
    active = active_r[]
    display(findall(active))
    if nps>=2 && method.verbose == true
        println("Identified np≥2 candidate phases. Moving on to stage III.")
        println("Candidate solutions are:")
        for i ∈ 1:nps
            xx = cache.ℳ[active]
            V = exp10.(first.(xx))
            x = getindex.(xx,Ref(2:nc+1))
            @show x,V
        end
    end
    if method.verbose == true
        println("=============================================")
        println("Stage III: Acceleration and convergence tests")
        println("=============================================")
        println("--------------------------------")
        println("Step 7: Free energy minimisation")
        println("--------------------------------")
    end

    @show length(cache.ℳ)
    return zeros(2,2),zeros(2,2),NaN
    ℳ = cache.ℳ
    X0 = Float64[]
    for i in active

    end
    X0 = vec(reshape(ℳˢ[:,1:nc-1],(1,nps*(nc-1))))

    X0 = append!(X0,ℳˢ[:,end])
    g(x) = Obj_HELD_tp_flash(model,p,T,n,x,nps)

    #Default options, feel free to change any of those
    options = OptimizationOptions(; x_abstol=0.0, x_reltol=0.0, x_norm=x->norm(x, Inf),
    g_abstol=1e-8, g_reltol=0.0, g_norm=x->norm(x, Inf),
    f_limit=-Inf, f_abstol=0.0, f_reltol=0.0,
    nm_tol=1e-8, maxiter=10000, show_trace=false)




    r = Solvers.optimize(g,X0,LineSearch(Newton()),options)
    if method.verbose==true
        println(r)
    end

    if method.verbose == true
        println("------------------------")
        println("Step 8: Convergence test")
        println("------------------------")
    end


    X = Solvers.x_sol(r)
    G = g(X)

    x = reshape(X[1:nps*(nc-1)],(nps,nc-1))
    x = Clapeyron.Fractions.FractionVector.(eachrow(x))
    V = exp10.(X[nps*(nc-1)+1:nps*nc])
    # ϕ = X[nps*nc+1:nps*(nc+1)]
    # λ = X[nps*(nc+1)+1:end]
    # if any(abs.(λ).<method.eps_μ) & method.verbose==true
    #     println("Mass balance could not be satisfied.")
    # end
    test_G = UBDⱽ-G
    μ = VT_chemical_potential.(model,V,T,x)/R̄/T
    test_μ = [abs((μ[j][i]-μ[j+1][i])/μ[j][i]) for i ∈ 1:nc for j ∈ 1:nps-1]

    @show test_μ
    @show test_G
    #test_μ = [abs((μ[j][i]-μ[j+1][i])/μ[j][i])<method.eps_μ for i ∈ 1:nc for j ∈ 1:nps-1]
    println(x)
    println(G)
    if (test_G >=method.eps_g) & all(<(method.eps_μ),test_μ)
        if method.verbose == true
            println("HELD has successfully converged to a solution. Terminating algorithm.")
        end
        return (collect(x),collect(ϕ.*x),G)
    else
        if method.verbose == true
            println("HELD has failed to converged to a solution. Terminating algorithm.")
        end
        return (collect(x),collect(ϕ.*x),G)
    end
end

function HELD_stage_II(cache,UBDⱽ,method,k)
    model, p, T, n = cache.model, cache.p, cache.T, cache.z
    ℳ = cache.ℳ
    λᴸ, λᵁ = cache.λᴸ, cache.λᵁ
    nc = length(n)
    G = cache.G
    if method.verbose == true
        println("-------------------------------------------------------")
        println("Step 3: Solve the outer problem (OPₓᵥ) at iteration k="*string(k))
        println("-------------------------------------------------------")
    end

    OPₓᵥ = Model(HiGHS.Optimizer)
    set_optimizer_attribute(OPₓᵥ, "log_to_console", false)
    set_optimizer_attribute(OPₓᵥ, "output_flag", false)
    @variable(OPₓᵥ, v)
    @variable(OPₓᵥ, λ[1:nc-1])
    @constraint(OPₓᵥ,v<=UBDⱽ)
    @constraint(OPₓᵥ,[i ∈ 1:length(G)],v<=G[i]+∑(λ.*(n[1:nc-1] .-ℳ[i][1:nc-1])))
    @constraint(OPₓᵥ,[i ∈ 1:nc-1],λᴸ[i]<=λ[i]<=λᵁ[i])
    @objective(OPₓᵥ, Max, v)
    optimize!(OPₓᵥ)
    λˢ = JuMP.value.(λ)
    UBDⱽ = JuMP.value.(v)
    if method.verbose == true
        println("UBDⱽ = "*string(UBDⱽ))
    end
    if method.verbose == true
        println("-------------------------------------------------------")
        println("Step 4: Solve the inner problem (IPₓᵥ) at iteration k="*string(k))
        println("-------------------------------------------------------")
    end
    i = 0
    while i < method.inner_steps
        Lⱽ(x) = _Lⱽ(model,p,T,x,λˢ,n)
        x0 = rand(length(n))
        x0 .= x0./sum(x0)
        v0 = log10(volume(model,p,T,x0))
        x0 = prepend!(x0,v0)
        r = Solvers.optimize(Lⱽ,x0[1:end-1])
        xᵏ = Solvers.x_sol(r)
        prepend!(xᵏ)
        Lⱽᵏ = Lⱽ(xᵏ)
        push!(xᵏ,1-sum(@view(xᵏ[2:end])))
        update_bounds!(cache,λˢ)
        if !already_in_cache(cache,xᵏ,true)
            add_candidate!(cache,xᵏ,precalc = true)
            cache.LV[end] = Lⱽᵏ
        end
        if Lⱽᵏ<UBDⱽ
            break
        end
        i+=1
    end
    if method.verbose == true
        println("Lⱽᵏ = "*string(cache.LV[end]))
        println("ℳ = "*string(ℳ[end]))
    end
    if method.verbose == true
        println("------------------------------------------------")
        println("Step 5: Select candidate phases at iteration k="*string(k))
        println("------------------------------------------------")
    end
    active = active_phases(cache,method,UBDⱽ,λˢ,k)
    return UBDⱽ,cache,active
end

function _pack_fraction(model,V,T,z)
    lb_v = lb_volume(model,z)
    ∑z = sum(z)
    return lb_v/V/sum(z)
end

function _pack_fraction(cache::HELDIPCache,i::Int)
    model,p,T = cache.model,cache.p,cache.T
    xi = cache.ℳ[i]
    V = first(xi)
    x = @view(xi[2:end])
    return _pack_fraction(model,V,T,x)
end

function active_phases(cache::HELDIPCache,method,UBDⱽ,λˢ,k)
    lenℳ = length(cache)

    #result. active phases to look
    active = zeros(Bool,lenℳ)
    #UBD test
    test_b = zeros(Bool,lenℳ)
    test_b2 = zeros(lenℳ)
    #λ test
    test_λ = zeros(Bool,lenℳ)
    test_λ2 = zeros(lenℳ)

    ℳ = cache.ℳ
    #η test
    test_η = ones(Bool,lenℳ)

    #fractions test
    test_x = ones(Bool,lenℳ)

    LV = cache.LV
    T = cache.T
    nc = length(cache.model)
    for i in 1:lenℳ

        #test that the difference between the current upper bound and bound of the phase is less than eps_b
        test_b[i] = abs(UBDⱽ- LV[i]) <= method.eps_b
        test_b2[i] = abs(UBDⱽ-LV[i])
        #test that the dual bonds are within tolerance
        μᵢ = cache.μ[i]
        Δλᵢ = abs.((μᵢ[1:nc-1] .- λˢ) ./ λˢ )
        test_λ[i] = norm(Δλᵢ,Inf) <= method.eps_λ
        test_λ2[i] = norm(Δλᵢ,Inf)

        xᵢ = @view(ℳ[i][2:end])
        ηᵢ = _pack_fraction(cache,i)
        for j in (i+1):lenℳ
            #we suppose that the last element has the the best LV
            xⱼ = @view(ℳ[i][2:end])
            ηⱼ = _pack_fraction(cache,j)

            #test that there are differences between each phase (volume)
            test_ηᵢⱼ = abs(ηᵢ-ηⱼ) >= method.eps_η
            test_η[i] = test_η[i] & test_ηᵢⱼ

            #test that there are differences between each phase (composition)
            test_xᵢⱼ = dnorm(xᵢ,xⱼ,Inf) >= method.eps_x
            test_x[i] = test_x[i] & test_xᵢⱼ
        end
    end
    if k == 500
        ii =  findall(test_b)
        @show ℳ[ii]
        @show findall(test_λ)
        @show findall(test_η)
        @show findall(test_x)
    end
    for i in 1:lenℳ
        active[i] = test_b[i] & test_λ[i] & (test_x[i] | test_η[i])
    end

    return active
end

function initial_candidate_phases(model,p,T,n,μ₀)
    nc = length(n)
    #=
   # x̂ = zeros(nc-1,nc+1)
    x̂ = [fill(0.0,nc) for i in 1:nc-1]
    #x̄ = zeros(nc-1,nc+1)
    x̄ = [fill(0.0,nc) for i in 1:nc-1]
    for i ∈ 1:nc-1
        x̂i = x̂[i]
        x̄i = x̄[i]
        x̂i[i] = n[i]/2
        x̄i[i] = (1+n[i])/2
        for k ∈ 1:nc-1
            if k != i
                x̂i[k] = (1-x̂i[i])/(nc-1)
                x̄i[k] = (1-x̄i[i])/(nc-1)
            end
        end
    end =#

    #x = vcat(x̂,x̄)
    #@show x
    x = Vector{Vector{Float64}}(undef,0)
    for i in 1:nc
        xi_pure = zeros(Float64,nc)
        xi_pure[i] = 1.
        push!(x,xi_pure)
    end

    #Fraction algebra
    for k in 2:4
        xk = Fractions.mul(n,k)
        xk_inverse = Fractions.mul(n,1/k)
        xkm = Fractions.neg(xk)
        xkm_inverse = Fractions.neg(n)
        push!(x,xk)
        push!(x,xkm)
        push!(x,xk_inverse)
        push!(x,xkm_inverse)
    end
    unique!(x)

    λ₀ = (μ₀[1:nc-1].-μ₀[nc])/R̄/T
    resize!(λ₀,nc-1)
    cache = HELDIPCache(model,p,T,n)
    for xi in x
        add_candidate_with_λ!(cache,xi,λ₀)
    end
    return cache
end

function Obj_HELD_tp_flash(model,p,T,x₀,X,np)
    nc = length(x₀)
    x = reshape(X[1:np*(nc-1)],(np,nc-1))

    V = exp10.(X[np*(nc-1)+1:np*nc])
    if np == 2
        ϕ = [(x₀[1]-x[1,1])/(x[2,1]-x[1,1]),(x₀[1]-x[2,1])/(x[1,1]-x[2,1])]
    elseif np==3
        ϕ = [0,
            ((x₀[1]-x[1,1])/(x[3,1]-x[1,1])-(x₀[2]-x[1,2])/(x[3,2]-x[1,2]))/((x[2,1]-x[1,1])/(x[3,1]-x[1,1])-(x[2,2]-x[1,2])/(x[3,2]-x[1,2])),
            ((x₀[1]-x[1,1])/(x[2,1]-x[1,1])-(x₀[2]-x[1,2])/(x[2,2]-x[1,2]))/((x[3,1]-x[1,1])/(x[2,1]-x[1,1])-(x[3,2]-x[1,2])/(x[2,2]-x[1,2]))]
        ϕ[1] = 1-sum(ϕ[2:3])
    end

    x = Clapeyron.Fractions.FractionVector.(eachrow(x))

    A = Clapeyron.eos.(model,V,T,x)
    f = sum(ϕ.*(A+p*V)./sum.(x))/Clapeyron.R̄/T
    F = f
    return F
end

export HELDTPFlash

#=
test_b = zeros(length(G))
    test_λ = zeros(length(G))
    test_cross_η = float.(LV.>LV')
    test_cross_x = zeros((length(G),length(G)))
    for m ∈ 1:length(G)
        test_b[m] += (UBDⱽ-LV[m]<=method.eps_b/R̄/T)
        vᵢ = exp10(ℳ[m,end])
        wᵢ = @view(ℳ[m,1:nc])
        μᵢ = VT_chemical_potential(model,vᵢ,T,wᵢ)/R̄/T
        μᵢ = μᵢ[1:nc-1].-μᵢ[nc]
        test_λ[m] += min(maximum(abs.((μᵢ[1:nc-1].-λˢ)./λˢ).>=method.eps_λ),1)
        ηm = packing_fraction(model,vᵢ,T,wᵢ)
        xm = ℳ[m,1:nc-1]
        for n ∈ 1:length(G)
            if n!=m
                vn = exp10(ℳ[n,end])
                wn = @view(ℳ[n,1:nc])
                ηn = packing_fraction(model,vn,T,wn)
                test_cross_η[m,n] += abs(ηm-ηn)<=method.eps_η
                test_cross_η[n,m] = test_cross_η[m,n]
                xn = ℳ[n,1:nc-1]
                test_cross_x[m,n] += dnorm(xm,xn,Inf)<=method.eps_x
                test_cross_x[n,m] = test_cross_x[m,n]
            end
        end
    end
    test = test_b+test_λ
    test_cross=test_cross_η+test_cross_x
    display(test_cross)
    test_cross = float.(test_cross.>=2)
    test .+= sum(test_cross,dims=2)
    test = Bool.(1 .-(test.>=1))
    display(test)
    ℳˢ = ℳ[test,:]
    Gˢ = G[test]
    LVˢ = LV[test]

=#