using JuMP, HiGHS
import Distributions
const HELD_verbose_dict = Dict{Int,String}(
    1001 => "Stability test and initialisation",
    1 => "Stability test at n₀",
    2 => "Initialisation of dual problem",
    1002 => "Solve the outer problem (OPx,V)"
)

_HELD_verbose_dict = Dict{Int,String}(
    1001 => "Stability test and initialisation",
    1 => "Stability test at n₀",
    2 => "Initialisation of dual problem",
    1002 => "Identification of candidate stable phases",
    3 => "Solve the outer problem (OPₓᵥ)",
    4 => "Solve the inner problem (IPₓᵥ)",
    5 => "Search for candidate stable phases",
    6 => "Increment iteration counter k = k + 1 and go to Step 3",
    1003 => "Acceleration and convergence tests",
    7 => "Minimisation of the Gibbs free energy over all candidate phases",
    8 => "Convergence test",
    9 => "Check for trace components.",
)


function held_verbose(stage,step,text = "")
    if step == 0
        data = _HELD_verbose_dict[stage + 1000]
        stage_str = repeat('I',stage)
        print(info_color("[ HELD - Stage " * stage_str * ": " * data))
    elseif stage == 0
        data = _HELD_verbose_dict[step]
        print(info_color("[ HELD - Step " * string(step) * ": " * data))
    end
    if text != ""
        print(info_color(" - "))
    end
    println(text)
end

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
# lb/V, t = 1, V = lb_V, t = 0 , V = inf
function _Lⱽ(model,p,T,x,λˢ,n,lb) #x = [eta,w]
    V = lb/x[1]
    w = Fractions.FractionVector(@view(x[2:end]))
    nc = length(model)
    (eos(model,V,T,w)+p*V)/R̄/T + sum(λˢ[j]*(n[j]-x[j+1]) for j ∈ 1:nc-1)
end


function _Lⱽ(model,p,T,x,λˢ,n) #x = [eta,w]
    w = Fractions.FractionVector(@view(x[2:end]))
    lb = lb_volume(model,w)
    V = lb/x[1]
    nc = length(model)
    (eos(model,V,T,w)+p*V)/R̄/T + sum(λˢ[j]*(n[j]-x[j+1]) for j ∈ 1:nc-1)
end

#Inner problem cache
struct HELDIPCache{𝕋}
    model::𝕋 #current eos model
    p::Float64 #pressure of the system, pa
    T::Float64 #temperature of the system, K
    z::Vector{Float64} #composition of the feed
    lb_v::Vector{Float64} #lower bound volumes
    ℳ::Vector{Vector{Float64}} #list of phases in the form [eta;w]
    G::Vector{Float64} #gibbs energies of the phases
    λᴸ::Vector{Float64} #lower bound of lagrangians
    λᵁ::Vector{Float64} #upper bound of lagrangians
    μ::Vector{Vector{Float64}} #list of chemical potentials
    LV::Vector{Float64} #list of tpd for each phase
    active::Vector{Bool} #list of active phases
end

Base.length(cache::HELDIPCache) = length(cache.G)

#constructor
function HELDIPCache(model,p,T,z)
    nc = length(z)
    lb_v = Vector{Float64}(undef,0)
    ℳ = Vector{Vector{Float64}}(undef,0)
    G = Vector{Float64}(undef,0)
    λᴸ = fill(Inf,nc - 1)
    λᵁ = fill(-Inf,nc - 1)
    μ = Vector{Vector{Float64}}(undef,0)
    LV = Vector{Float64}(undef,0)
    active = Vector{Bool}(undef,0)
    return HELDIPCache(model,Float64(p),Float64(T),z,lb_v,ℳ,G,λᴸ,λᵁ,μ,LV,active)
end

function add_candidate!(cache::HELDIPCache,xx;precalc = false)
    p,T = cache.p,cache.T
    model = cache.model
    nc = length(model)
    if precalc
        x = xx[2:end]
        @show xx
        lb_v_new = lb_volume(model,x)
        η = first(x)
        #η = lb_v/v
        v = lb_v_new/η
        m = xx
    else
        @show xx
        x = xx
        lb_v_new = lb_volume(model,x)
        v = volume(model,p,T,xx)
        η = lb_v_new/v
        m = vcat(η,x)
    end
    μ_new = VT_chemical_potential(model,v,T,x)/R̄/T

    G_new = dot(μ_new,x)
    if isnan(G_new)
        G_new =  VT_gibbs_free_energy(model,v,T,x)/R̄/T
    end
    μλ = copy(μ_new)
    μλ_nc = μλ[end]
    μλ .-= μλ_nc
    resize!(μλ,nc-1)
    nc = length(cache.model)
    #adding the calculated values:
    push!(cache.lb_v,lb_v_new)
    push!(cache.ℳ,m)
    push!(cache.G,G_new)
    update_bounds!(cache::HELDIPCache,μλ)
    push!(cache.μ,μ_new)
    push!(cache.LV,1.)

    #by default a new phase is active. after setting the phase to not active
    #the decision is final, but it needs to be calculated in another function.
    push!(cache.active,true)
    return cache
end

function already_in_cache(cache::HELDIPCache,xx,method,precalc = false)
    model = cache.model
    T = cache.T
    p = cache.p
    if !precalc
        v = volume(model,p,T,xx)
        x = xx
        lb = lb_volume(model,x)
        m = vcat(lb/v,x)
    else
        m = xx
    end
    x₀ = @view(m[2:end])
    @show x₀
    for xi in cache.ℳ
        xᵢ = @view(xi[2:end])
        @show xᵢ
        if dnorm(x₀,xᵢ,Inf) <= 1e-8
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

function update_bounds!(cache::HELDIPCache)
    nc = length(cache.model)
    for μi in cache.μ
        μλ = copy(μi)
        μλ_nc = μλ[end]
        μλ .-= μλ_nc
        resize!(μλ,nc-1)
        update_bounds!(cache,μλ)
    end
end

function add_candidate_with_λ!(cache::HELDIPCache,x,λ)
    #@show λ
    cache = add_candidate!(cache,x)
    xx = cache.ℳ[end]
    model,p,T = cache.model,cache.p,cache.T
    n = cache.z
    lb = cache.lb_v[end]
    cache.LV[end] = 1
    return cache,true
end

function HELD_stage_I(model,p,T,n,method,cache,μ₀,v₀)
    nc = length(model)
    p̃ = method.p_tun
    TYPE = typeof(p+T+first(n))
    # Stage I: Stability test and initialisation
    if method.verbose == true
        held_verbose(1,0) #stage 1
        held_verbose(0,1) #step 1
    end
    #v00 = volume(model,p,T,n)
    vₛ = log10(v₀)
    μ₀ = VT_chemical_potential(model,v₀,T,n)
    xₛ = prepend!(deepcopy(n),vₛ)

    function d(x)
        vv = exp10(x[1])
        xx = @view(x[2:end])
        return (eos(model,vv,T,xx)+p*vv)/sum(xx)-dot(xx,μ₀)
    end
    method.verbose == true && held_verbose(0,1,"Beginning tunneling")

    i = 0
    while i<10nc
        function f(x)
            dx = d(x)
            xx = @view(x[2:end])
            xxₛ = @view(xₛ[2:end])
            return dx*exp(p̃/dnorm(xx,xxₛ))
        end
        #f(x) = d(x)*exp(p̃/√(∑((xₛ[2:end]-x[2:end]).^2)))
        x0 = rand(length(n))
        x0 = x0./sum(x0)
        v0 = log10(volume(model,p,T,x0))
        x0 = prepend!(x0,v0)
        r = Solvers.optimize(f,x0)
        xₙ = Solvers.x_sol(r)
        dₙ = d(xₙ)
        if dₙ<0
            if method.verbose == true
                held_verbose(0,1,"Negative tangent found. Moving on to step 2.")
            end
            xnew = xₙ[2:end] 
            xnew ./= sum(xnew)
            add_candidate!(cache,xnew,precalc = false)
            cache.LV[end] = 1
            return true
        end
        i+=1
    end
    if i==10nc
        if method.verbose == true
            held_verbose(0,1,"No negative tangent found. Initial point is stable.")
            held_verbose(0,1,error_color("Terminating HELD"))
        end
        return false
    end

    if method.verbose == true
        held_verbose(0,2)
    end
    
end

function tp_flash_impl(model::EoSModel, p, T, n, method::HELDTPFlash)
    UBDⱽ = gibbs_free_energy(model,p,T,n)/R̄/T
    v00 = volume(model,p,T,n)
    μ₀ = VT_chemical_potential(model,v00,T,n)
    cache = initial_candidate_phases(model,p,T,n,μ₀)
    verbose = method.verbose
    #λ₀ = (μ₀[1:nc-1].-μ₀[nc])/R̄/T
    #LV = vec(G.+ sum(λ₀'.*(n[1:nc-1]'.-ℳ[:,1:nc-1]),dims=2))
    tunnelling_success = HELD_stage_I(model,p,T,n,method,cache,μ₀,v00)
    !tunnelling_success && return cache
    k = 0
    if verbose == true
        held_verbose(0,2,"Iteration counter set to k="*string(k))
        held_verbose(0,2,"Upper bound set to UBDⱽ="*string(UBDⱽ))
        held_verbose(2,0)
    end
    nps=1
    k = 0
    valid_phs = 0
    proceed_with_next_step = false
    suggestion = :reuse_starting
    while k<=method.outer_steps
        UBDⱽ,cache,λˢ = HELD_stage_II(cache,UBDⱽ,method,k,suggestion)
        suggestion = :none #reset suggestion
        active = active_phases!(cache,method,UBDⱽ,λˢ,k)
        clean_phases!(cache,method)
        if count(active) >=2
            proceed_with_next_step = true
        end

        #if count(active) == 1
        #    #suggestion = :complete_simplex
        #end

        #step 7: optimization of gibbs energy
        if proceed_with_next_step
            if verbose
                held_verbose(3,0)
                held_verbose(0,7)
            end
            #preprocessing: see if current feed inside the simplex formed by the candidate phases
            #clean_phases!(cache,method)
            #in_simplex = feed_in_simplex(cache,verbose)
            
            #if !in_simplex
            #    verbose && held_verbose(0,7,"Candidate phases don't from a simplex that contains the feed")
            #    verbose && held_verbose(0,7,"trying to clean phases")
            #    in_simplex = clean_phases!(cache,method)
            #end

            #if !in_simplex
            #    verbose && held_verbose(0,7,"cleaning unsucessful. Returning to stage II.")
            #    proceed_with_next_step = false
            #else
                #verbose && held_verbose(0,7,"Candidate phases from a sucessful simplex. proceeding with Gibbs optimization")
            
            
            sucess,cache =  HELD_stage_III(cache,method)
            proceed_with_next_step = sucess
            @show sucess
            return cache
            # end
        end

        #step 8: convergence test
        if proceed_with_next_step
            converged = HELD_convergence_test(cache,method,UBDⱽ)
            if !converged
                proceed_with_next_step = false
            else
                #we found a current minimum
                #TODO: "obtain the true eq for any trace comps"
                break
            end
        end

        valid_phs = 0
        k+=1
    end

    return cache
end

function HELD_stage_II(cache,UBDⱽ,method,k,suggestion = :none)
    model, p, T, n = cache.model, cache.p, cache.T, cache.z
    ℳ = cache.ℳ
    λᴸ, λᵁ = cache.λᴸ, cache.λᵁ
    nc = length(n)
    G = cache.G
    if method.verbose == true
        held_verbose(0,3,"iteration $k")
    end

    OPₓᵥ = Model(HiGHS.Optimizer)
    set_optimizer_attribute(OPₓᵥ, "log_to_console", false)
    set_optimizer_attribute(OPₓᵥ, "output_flag", false)
    
    @variable(OPₓᵥ, v)
    @variable(OPₓᵥ, λ[1:nc-1])
    
    @constraint(OPₓᵥ,v <= UBDⱽ)

    #v <= Gi + ∑(λi*(n-xi))
    @constraint(OPₓᵥ,[i ∈ 1:length(G)],v <= G[i]+sum(λ.*(n[1:nc-1] .- ℳ[i][2:nc-1])))
    
    #λᴸ <= λ <= λᵁ 
    @constraint(OPₓᵥ,[i ∈ 1:nc-1],λᴸ[i] <= λ[i] <= λᵁ[i])
    @objective(OPₓᵥ, Max, v)
    optimize!(OPₓᵥ)
    λˢ = JuMP.value.(λ)
    new_UBDⱽ = JuMP.value.(v)
    if method.verbose == true
        held_verbose(0,3,"UBDⱽ = "*string(UBDⱽ))
        #println("UBDⱽ = "*string(UBDⱽ))
    end
    if method.verbose == true
        held_verbose(0,4)
    end
    
    if suggestion == :reuse_starting
        max_inner = length(ℳ)
    else
        max_inner = method.inner_steps
    end
    i = 0
    while i < max_inner
        Lⱽ(x) = _Lⱽ(model,p,T,x,λˢ,n)


        #if suggestion == :complete_simplex
            #x0 = complete_simplex_suggestion(cache)
        if suggestion == :reuse_starting #
            x0 = ℳ[i+1]
        else
            x_new = rand(length(n))
            x_new .= x_new./sum(x_new)
            lb = lb_volume(model,x_new)
            ln_lb = log(lb)
            rn = ln_lb + rand()*(1 - ln_lb)
            credible_eta = min(1,lb/exp(rn))
            x_new = prepend!(x_new,credible_eta)
            x0 = x_new
        end
        
        r = Solvers.optimize(Lⱽ,x0[1:end-1])
        xᵏ = Solvers.x_sol(r)
        prepend!(xᵏ)
        Lⱽᵏ = Lⱽ(xᵏ)
        push!(xᵏ,1-sum(@view(xᵏ[2:end])))
        update_bounds!(cache,λˢ)
        #if !already_in_cache(cache,xᵏ,method,true)
        #if suggestion == :complete_simplex
        #    if !already_in_cache(cache,xᵏ,method,true)
                add_candidate!(cache,xᵏ,precalc = true)
        #    end
        #end
        #end
        #@show Lⱽᵏ
        if Lⱽᵏ < new_UBDⱽ
            if !already_in_cache(cache,xᵏ,method,true)
                held_verbose(0,4,"we add a new xᵏ = $xᵏ")
                add_candidate!(cache,xᵏ,precalc = true)
                cache.LV[end] = 1
            end
            new_UBDⱽ = Lⱽᵏ
            if suggestion != :reuse_starting
        
                break
            end
            #return Lⱽᵏ,cache,λˢ
        end
        #if suggestion == :complete_simplex
        #    suggestion = :none
        #end
        i+=1
    end

    if method.verbose == true
        held_verbose(0,5,"iteration $k")
    end

    return new_UBDⱽ,cache,λˢ
end

function _pack_fraction(model,V,T,z)
    lb_v = lb_volume(model,z)
    ∑z = sum(z)
    return lb_v/V/sum(z)
end

function _pack_fraction(cache::HELDIPCache,i::Int)
    lb_v = cache.lb_v[i]
    η = cache.ℳ[i][1]
    return l
end

function active_phases!(cache::HELDIPCache,method,UBDⱽ,λˢ,k)
    active = cache.active
    ℳ = cache.ℳ
    μ = cache.μ
    nc = length(cache.model)
    n = cache.z
    G = cache.G
    for i in 1:length(active)
        xᵢ = @view(ℳ[i][2:end])
        wᵢ = Fractions.FractionVector(xᵢ)
        LVˢᵢ = G[i] + sum(λˢ[j]*(n[j]-wᵢ[j+1]) for j ∈ 1:nc-1)
        test_b = UBDⱽ - LVˢᵢ <= method.eps_b
        active[i] = test_b
        μᵢ = μ[i]
        Δλ = -Inf
        for k in 1:nc-1
            Δλ = max(Δλ,abs((μᵢ[k] -  λˢ[k])/λˢ[k]))
        end
        test_λ = Δλ <= method.eps_λ
        test_initial = cache.LV[i] < 0
        active[i] = test_b & test_λ & test_initial
        #we perform the checks on each active phase.
        #not active phases are discarded
        if active[i]
            #test that the difference between the current upper bound and bound of the phase is less than eps_b
            ηᵢ = ℳ[i][1]
            for j in 1:i-1
                if active[j]
                    xⱼ = @view(ℳ[i][2:end])
                    ηⱼ = ℳ[j][1]
                    test_x = dnorm(xᵢ,xⱼ,Inf) >= method.eps_x
                    test_η = abs(ηᵢ-ηⱼ) >= method.eps_η

                    #if the phases are the same, deactivate one
                    if !test_x & !test_η
                        active[j] = false
                    end
                end
            end
        end
    end
    return active
end

#mixes all the existing phases into one, creates a second phase that by defintion completes the simplex.
function complete_simplex_suggestion(cache::HELDIPCache)
    active = cache.active
    ℳ = cache.ℳ
    active_phases = ℳ[active]
    n = cache.z
    nc = length(cache.model)
    #we obtain a "mean phase" consisting in a mixing of all existing phases
    w = zeros(Float64,nc + 1)
    for xi in active_phases
        w .+= xi
    end
    ∑mean = sum(w)
    w ./= ∑mean
    w = @view(w[2:end])
    
    nmin,nmax = extrema(n)
    wmin,wmax = extrema(n)
    #=
    zi - wi * β  = 0
    minimum case:
    zmax - wmin*β = 0 => β = zmax/wmin
    maximum case
    zmin - wmax*β = 0 => β = zmin/wmax
    =#
    β1 =  nmax/wmin
    β2 =  nmin/wmax
    if 0 < β1 < 1
        β = β1
    else
        β = β2
    end
    #     wmin
    y = (n .- β .* w) ./(1 .- β)
    y ./= sum(y)

    model = cache.model
    T = cache.T
    p = cache.p
    lb_v = lb_volume(model,y)
    v = volume(model,p,T,y)
    prepend!(y,lb_v/v)
    return y
end

function initial_candidate_phases(model,p,T,n,μ₀)
    nc = length(n)
    x = initial_candidate_fractions(n)
    unique!(x)
    λ₀ = (μ₀[1:nc-1].-μ₀[nc])/R̄/T
    resize!(λ₀,nc-1)
    cache = HELDIPCache(model,p,T,n)
    for xi in x
        add_candidate_with_λ!(cache,xi,λ₀)
    end
    cache.LV .= -1
    return cache
end



function Obj_HELD_tp_flash(model,p,T,x₀,X,lb)
    np = length(lb)
    nc = length(x₀)
    x = reshape(@view(X[1:np*(nc-1)]),(nc-1,np))

    #x = reshape(@view(X[1:np*(nc-1)]),(np,nc-1))
    #x = xx
    #display(x)
    V = lb ./ (X[(np*(nc-1)+1):end])
    #=
    if np == 2
        ϕ = [(x₀[1]-x[1,1])/(x[2,1]-x[1,1]),(x₀[1]-x[2,1])/(x[1,1]-x[2,1])]
    elseif np==3
        ϕ = [0,
            ((x₀[1]-x[1,1])/(x[3,1]-x[1,1])-(x₀[2]-x[1,2])/(x[3,2]-x[1,2]))/((x[2,1]-x[1,1])/(x[3,1]-x[1,1])-(x[2,2]-x[1,2])/(x[3,2]-x[1,2])),
            ((x₀[1]-x[1,1])/(x[2,1]-x[1,1])-(x₀[2]-x[1,2])/(x[2,2]-x[1,2]))/((x[3,1]-x[1,1])/(x[2,1]-x[1,1])-(x[3,2]-x[1,2])/(x[2,2]-x[1,2]))]
        ϕ[1] = 1-sum(ϕ[2:3])
    end    
    @show ϕ =#
    z0 = zeros(eltype(X),nc)
    v0 = zeros(eltype(X),nc)
    A = zeros(eltype(X),(nc,np - 1))
    for (i,xi) in pairs(eachcol(x))
        if i == 1
            w0  = FractionVector(xi)
            v0 .= w0
            z0 .= x₀ .- v0
        else
            ai = @view(A[:,i-1])
            vi  = FractionVector(xi)
            ai .= vi .- v0
        end
    end
    try
        ϕ = A\z0
        prepend!(ϕ,1-sum(ϕ))
        f = zero(typeof(p+T+first(x₀)))
        for (i,zi) in pairs(eachcol(x))
            xi = FractionVector(zi)
            Ai = eos(model,V[i],T,xi)
            f += ϕ[i]*(Ai + p*V[i])/sum(xi)
        end
        f = f/R̄/T
        return f
    catch
        _0 = zero(eltype(x))
        return _0/_0
    end
end

export HELDTPFlash

function feed_in_simplex(cache::HELDIPCache,verbose = false)
    α = feed_simplex(cache,cache.z)
    α0 = 1-sum(α)
    α = vcat(α0,α)
    verbose && held_verbose(0,7,"simplex result: $α")
  
    #check simplex conditions. #TODO: allow for aproximate search?
    in_simplex = (sum(α) <= 1) & all(>(0),α)
    return in_simplex
end

function feed_simplex(cache::HELDIPCache,x0 = cache.z)
    nc = length(x0)
    active = cache.active
    ℳ = cache.ℳ

    #we need to check if the candidate compositions (or a subset of those)
    #this is equivalent to say: x0 is inside a simplex made by the candidate solutions
    possible_sols = ℳ[active]
    l = length(possible_sols)
    #check for l-simplex
    A = zeros(nc,l - 1) #pivoted values
    z0 = zeros(nc) #pivoted feed
    i0 = 0 #pivot
    v0 = @view(possible_sols[1][2:end])
    z0 .= x0 .- v0
    #build pivoted matrix
    for i in 1:l-1
        ai = @view(A[:,i])
        vi = @view(possible_sols[i+1][2:end])
        ai .= vi .- v0
    end
    #solve the system of equations
    try
        α = A\z0
        return α
    catch
        α = fill(NaN,l - 1)
        return α
    end
end

function initial_candidate_fractions(n)
    nc = length(n)
    x̂ = [zeros(nc) for i in 1:nc-1]
    x̄ = [zeros(nc) for i in 1:nc-1]
    xp = Vector{Vector{Float64}}(undef,0)
    
    #generation by the method in Pereira et al. (2010).
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
        x̂i[nc] = 1 - sum(x̂i)
        x̄i[nc] = 1 - sum(x̄i)
    end

    #pure generation
    k = 100000
    for i in 1:nc
        xi = fill(1.,nc)
        xi[i] = k
        sumxi = sum(xi)
        xi ./= sumxi
        push!(xp,xi)
        for j in (0.1,0.5,0.7)
            push!(xp,Fractions.mul(xi,j))
        end
    end
    x = vcat(x̂,x̄)
    return x
end

#free gibbs energy minimization
function HELD_stage_III(cache::HELDIPCache,method)
    active = cache.active
    ℳˢ = cache.ℳ[cache.active]
    i1 = 1:length(cache.active)
    nc = length(cache.model)
    np = min(nc,4) #max of 4 simultaneous phases.
    idx = sortperm(cache.G[active],rev = true) #sort by amount of gibbs energy
    idx1 = i1[idx][1:np]

    α =  feed_simplex(cache)
    #prepend!(α,1-sum(α))
    ℳˢ̄ = ℳˢ[idx][1:np] #the phases that have the least amount of gibbs energy.    
    model,p,T,n = cache.model,cache.p,cache.T,cache.z
    
    #build initial vectors
    X0 = reduce(vcat,(@view(mi[2:end-1]) for mi in ℳˢ̄))
    V0 = [first(mi) for mi in ℳˢ̄]
    lb = cache.lb_v[cache.active][idx][1:np]
    append!(X0,V0)
    


    #X0 = append!(X0,ℳˢ[:,end])
    g(x) = Obj_HELD_tp_flash(model,p,T,n,x,lb)

    #Default options, feel free to change any of those
    options = OptimizationOptions(; x_abstol=0.0, x_reltol=0.0, x_norm=x->norm(x, Inf),
    g_abstol=1e-8, g_reltol=0.0, g_norm=x->norm(x, Inf),
    f_limit=-Inf, f_abstol=0.0, f_reltol=0.0,
    nm_tol=1e-8, maxiter=100, show_trace=false)

    r = Solvers.optimize(g,X0,LineSearch(Newton()),options)
    #update result with the current values
    r_min = Solvers.x_minimum(r)
    x_r = Solvers.x_sol(r)
    @show x_r
    if isnan(r_min)
        #failed to converge, #TODO: do something about that.
        #maybe merge the phases? delete those?
        println(error_color("fail in Newton Optimization. merging phases"))
        return false,cache
    else
        #@show active
        #active .= false
        #active[idx1] .= true
    
        x = reshape(@view(x_r[1:np*(nc-1)]),(nc-1,np))
        @show collect(eachcol(x))

        if sum(x)/length(x)/first(x) ≈ 1.
            throw(error("all are the same!"))
        end
        #x = reshape(@view(X[1:np*(nc-1)]),(np,nc-1))
        #x = xx
        #display(x)
        eta = x_r[(np*(nc-1)+1):end]
        #update phases to results of optimization.
        ki = 0
        for i in 1:np
            xi = @view(x[:,1])
            wi = vcat(xi,1 - sum(xi))
            if !already_in_cache(cache,wi,method,false)
                Vi = lb[i]/eta[i] #modified eta
                lb_i = lb_volume(model,wi)
                eta_i = lb_i/Vi
                prepend!(wi,eta_i) 
                add_candidate!(cache,wi;precalc = true)
            end
        end
        update_bounds!(cache::HELDIPCache)
        return true,cache
    end
end

function HELD_convergence_test(cache::HELDIPCache,method,UBDⱽ)
    active = cache.active
    ℳ = cache.ℳ
    μ = cache.μ
    nc = length(cache.model)
    n = cache.z
    G = cache.G
    convergence_test = false
    test_G = true
    test_μ = true
    i_0 = 0 #first 
    for i in 1:length(active)
         #we perform the checks on each active phase.
        #not active phases are discarded
        if active[i]
            #test that G is less than, but close to the UBDⱽ
            test_G = test_G & (0 <= G[i] - UBDⱽ <= method.eps_g)
            μi = μ[i]
            j = i+1
            while j <= length(active)
                if active[j]
                    μj = μ[j]
                    for k in 1:nc
                        test_μ = test_μ & (abs((μi[k] - μj[k])/μi[k]) < method.eps_μ)  
                    end
                    break
                end
                j += 1
            end
        end
    end
    convergence_test = test_G & test_μ
    return convergence_test
end

function Base.deleteat!(cache::HELDIPCache,v)
    deleteat!(cache.lb_v,v)
    deleteat!(cache.ℳ,v)
    deleteat!(cache.G,v)
    deleteat!(cache.μ,v)
    deleteat!(cache.LV,v)
    deleteat!(cache.active,v)
end

function clean_phases!(cache::HELDIPCache,method)
    iall = 1:length(cache.active)
    #alpha = feed_simplex(cache)
    delete_idx = Int[]
    for i in iall
        mi = @view cache.ℳ[i][2:end]
        for j in (i+1):length(cache.active)
            mj = @view cache.ℳ[j][2:end]
            if dnorm(mi,mj,Inf) < method.eps_x
                push!(delete_idx,j)
            end
        end
    end
    sort!(unique!(delete_idx))
    @show delete_idx
    deleteat!(cache,delete_idx)
    
    #alpha = feed_simplex(cache)
    #@show alpha
    #return feed_in_simplex(cache)
    return nothing
end
