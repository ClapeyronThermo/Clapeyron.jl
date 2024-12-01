struct Vfrac
    k::Int
end

struct FlashSpecifications{S1,V1,S2,V2}
    spec1::S1
    val1::V1
    spec2::S2
    val2::V2

    function FlashSpecifications(spec1::S1,val1::V1,spec2::S2,val2::V2) where {S1,V1,S2,V2}
        @assert spec1 == volume || spec1 == temperature || spec1 == entropy || spec1 == enthalpy || spec1 == internal_energy || spec1 == pressure || spec1 isa Vfrac
        @assert spec2 == volume || spec2 == temperature || spec2 == entropy || spec2 == enthalpy || spec2 == internal_energy || spec2 == pressure || spec2 isa Vfrac
        @assert spec1 !== spec2
        return new{S1,V1,S2,V2}(spec1,val1,spec2,val2)
    end
end

const FLASH_SPECS_FN = (
    volume,
    pressure,
    temperature,
    enthalpy,
    entropy,
    internal_energy,
    Vfrac(2),
)

function FlashSpecifications(;v = nothing,T = nothing,p = nothing,h = nothing,s = nothing,u = nothing,q = nothing)
    values = (v,p,T,h,s,u,q)
    nothing_values = isnothing.(values) .|> !
    if count(nothing_values) != 2
        throw(ArgumentError("invalid FlashSpecifications input argument, you need to specify exactly two specifications"))
    end
    i1,i2 = findfirst(nothing_values),findlast(nothing_values)
    val1,val2 = values[i1],values[i2]
    spec1,spec2 = FLASH_SPECS_FN[i1],FLASH_SPECS_FN[i2]
    return FlashSpecifications(spec1,val1,spec2,val2)
end

export FlashSpecifications

function Base.show(io::IO,x::FlashSpecifications)
    spec1,spec2,val1,val2 = x.spec1,x.spec2,x.val1,x.val2
    print(io,"FlashSpecifications(")
    print(io,spec1," => ",val1,", ")
    print(io,spec2," => ",val2,")")
end

spec_intensive(::typeof(pressure)) = true
spec_intensive(::typeof(temperature)) = true
spec_intensive(x) = false


function normalize_spec(s::FlashSpecifications,k)
    _1 = oneunit(1/k)
    spec1,spec2,val1,val2 = s.spec1,s.spec2,s.val1,s.val2
    if spec_intensive(spec1)
        newval1 = val1*_1
    else
        newval1 = val1/k
    end

    if spec_intensive(spec2)
        newval2 = val2*_1
    else
        newval2 = val2/k
    end

    return FlashSpecifications(spec1,newval1,spec2,newval2)
end

function set_vfrac(s::FlashSpecifications,i)
    if s.spec1 isa Vfrac
        return FlashSpecifications(Vfrac(i[s.spec1.k]),s.val1,s.spec2,s.val2)
    elseif s.spec2 isa Vfrac
        return FlashSpecifications(s.spec1,s.val1,Vfrac(s.spec2.k),s.val2)
    else
        return s
    end
end

function xy_input_to_flash_vars(input,np,nc,z)
    idx_comps_end = nc*np
    idx_comps = 1:idx_comps_end
    idx_volumes = (1:np) .+ idx_comps_end
    idx_β = (1:np) .+ idx_volumes[end]
    comps = @view input[idx_comps]
    volumes = @view input[idx_volumes]
    β = @view input[idx_β]
    #fill last component vector
    return comps,β,volumes
end

function xy_input_to_result(input,np,nc,z)
    compsx,βx,volumes = xy_input_to_flash_vars(input,np,nc,z)
    TT = eltype(input)
    comps = Vector{Vector{TT}}(undef,np)
    β = Vector{TT}(undef,np)
    for j in 1:np
        compsi = viewn(compsx,nc,j)
        comps[j] = collect(compsi)
        β[j] = βx[j]
    end
    T = input[end]
    #we return in order of increasing molar volumes
    idx = sortperm(volumes)
    return comps[idx],β[idx],volumes[idx],T
end

function spec_to_vt end

for prop in [:entropy,:enthalpy,:temperature,:pressure,:internal_energy,:gibbs_free_energy,:helmholtz_free_energy]
    @eval begin
        function spec_to_vt(model,V,T,z,spec::typeof($prop))
            VT.$prop(model,V,T,z)
        end
    end
end

spec_to_vt(model,V,T,z,spec::typeof(volume)) = V


#s = dadt
requires_pv(::typeof(entropy)) = false
requires_st(::typeof(entropy)) = true
requires_a(::typeof(entropy)) = false

#a = a
requires_pv(::typeof(helmholtz_free_energy)) = false
requires_st(::typeof(helmholtz_free_energy)) = false
requires_a(::typeof(helmholtz_free_energy)) = true

#g = a + pv
requires_pv(::typeof(gibbs_free_energy)) = true
requires_st(::typeof(gibbs_free_energy)) = false
requires_a(::typeof(gibbs_free_energy)) = true

#u = a + st
requires_pv(::typeof(internal_energy)) = false
requires_st(::typeof(internal_energy)) = true
requires_a(::typeof(internal_energy)) = true

#h = a + st + pv
requires_pv(::typeof(enthalpy)) = true
requires_st(::typeof(enthalpy)) = true
requires_a(::typeof(enthalpy)) = true

#v = v
requires_pv(::typeof(volume)) = true
requires_st(::typeof(volume)) = false
requires_a(::typeof(volume)) = false

#fallback: temperature, pressure, (fractions when added)
requires_pv(x) = false
requires_st(x) = false
requires_a(x) = false

requires_pv(x::FlashSpecifications) = requires_pv(x.spec1) | requires_pv(x.spec2)
requires_st(x::FlashSpecifications) = requires_st(x.spec1) | requires_st(x.spec2)
requires_a(x::FlashSpecifications) = requires_a(x.spec1) | requires_a(x.spec2)

struct XYFlashTag end

function xy_flash_neq(output,model,zbulk,np,input,state::F,μconfig) where F
    #=

    inspired in teqp and Ben Gharbia (2016)
    variables:
    phase fractions (β) [np]
    extended phase compositions (ξi) [np*nc]
    phase volumes (vi) [np] {note 1}
    T [1]
    total: np*nc + 2*np + 1
    ----------------------------------------
    equations:
    pressure equalities: [np - 1]
    chemical potential equalities: [(np-1)*nc]
    specification 1: [1]
    specification 2: [1]

    total: np*nc + np - nc + 1
    ----------------------------------------
    missing: np + nc
    ----------------------------------------
    restrictions:
    sum(βj*ξij) - zi = 0 [nc] {note 2}
    min(βj,1-sum(ξj)) = 0 [np] {note 3}

    ----------------------------------------
    domain:
    0 <= βj
    0 <= ξj
    ----------------------------------------
    specifications:
    spec = hset - sum(βi*hi)
    spec = sset - sum(βi*si)

    for caloric properties
    spec = (x0 - δa*∑βiai - δpv*p*∑βivi - δst*T*∑βisi)/RT
    δx := 1 if the term is required, 0 if its not.

    for pressure:
    spec = (p - p0)/p_scale

    for temperature (slack):
    spec = (T - T0)/T0

    for volume:
    spec = (v - ∑βivi)v_scale
    v_scale = RT/p_scale

    ----------------------------------------
    notes:
    {1} Ben Gharbia (2016) supposes that volumes are calculated each step.
        we iterate with the volumes. this adds np new equations, removes np volume calculations.
    {2} The first relation can be used to reduce the number of variables. the second one needs to be added as additional equations.
    {3} This relation allows simultaneous phase stability and phase equilibria calculation.

    =#
    T = input[end]
    nc = length(model)

    flash_vars = xy_input_to_flash_vars(input,np,nc,zbulk)
    ξ,β,volumes = flash_vars
    #@show primalval(ξ)
    #@show primalval(β)
    #@show primalval(volumes)
    β_end,j_end = β[np],np
    w_end = viewn(ξ,nc,j_end)
    v_end = volumes[j_end]
    needs_a = requires_a(state)
    needs_st = requires_st(state)
    needs_pv = requires_pv(state)
    
    
    if needs_st
        a_end,dadv_end,dadT_end = ∂f_vec(model,v_end,T,w_end)
        p_end,s_end = -dadv_end,-dadT_end
    elseif needs_a && !needs_st
        a_end,dadv_end = f∂fdV(model,v_end,T,w_end)
        p_end,s_end = -dadv_end,zero(a_end)
    else
        p_end = pressure(model,v_end,T,w_end)
        s_end,a_end = zero(p_end),zero(p_end)
    end
    #@show primalval(p_end)
    #@show primalval(pressure(model,v_end,T,w_end))
    #@show primalval(w_end)
    ∑a = β_end*a_end
    ∑s = β_end*s_end
    ∑v = needs_pv ? β_end*v_end : zero(∑a)
    
    ps = p_scale(model,zbulk)
    R = Rgas(model)
    RT = R*T
    RTinv = 1/RT
    vs = RT/ps

    #fill pressure constraints:
    idx_p_constraints = 1:(np-1)
    p_constraints = @view output[idx_p_constraints]
    jj = 0
    for j in 1:np
        j == j_end && continue
        jj += 1
        βj = β[j]
        vj = volumes[j]
        if needs_pv
            ∑v += βj*vj
        end
        wj = viewn(ξ,nc,jj)
        if needs_st
            aj,dadvj,dadTj = ∂f_vec(model,volumes[j],T,wj)
            pj,sj = -dadvj,-dadTj
        elseif needs_a && !needs_st
            aj,dadvj = f∂fdV(model,volumes[j],T,wj)
            pj,sj = -dadvj,zero(aj)
        else
            pj = pressure(model,volumes[j],T,wj)
            sj,aj = zero(pj),zero(pj)
        end
        #@show primalval(pj)
        ∑a += βj*aj
        ∑s += βj*sj
        p_constraints[jj] = (pj - p_end)/ps
    end
    #fill chemical potential equalities:

    idx_μ_constraints = (1:((np-1)*nc)) .+ (np - 1)
    μ_constraints = @view output[idx_μ_constraints]
    μ_end = similar(output,nc)

    VT_chemical_potential_res!(μ_end,model,v_end,T,w_end,μconfig)

    jj = 0
    for j in 1:np
        j == j_end && continue
        jj += 1
        Fj = @inbounds viewn(μ_constraints,nc,jj)
        for i in 1:nc
            Fj[i] = exp(μ_end[i]*RTinv)*w_end[i]/v_end
        end
    end
    jj = 0
    for j in 1:np
        j == j_end && continue
        jj += 1
        vj = @inbounds volumes[j]
        wj = viewn(ξ,nc,j)
        VT_chemical_potential_res!(μ_end,model,vj,T,wj)
        Fj = viewn(μ_constraints,nc,jj)

        for i in 1:nc
            μji = μ_end[i]
            #=
            Δuᵣ = μ1i - μji
            Δu = Δuᵣ*RTinv + log(vj*w1[i]/(v1*wj[i]))
            Fj[i] = Δu #exponentiating
            expΔu = vj*w1[i]*/(v1*wj[i])exp(μ1i*RTinv)/exp(μji*RTInv)
            1 = vj*w1[i]*/(v1*wj[i])exp(μ1i*RTinv)/exp(μji*RTInv)
            exp(μji*RTInv)*wj[i]/vj = exp(μ1i*RTinv)/exp(μji*RTInv)*w1[i]/v1
            =# 
            
            Fj[i] -= exp(μji*RTinv)*wj[i]/vj
            Fj[i] *= vs
        end
    end

    #fill β,extended composition constraints
    βξspec_constraints = viewlast(output,np+nc+2)
    ξ_constraints = @view βξspec_constraints[np+1:end-2]
    β_constraints = @view βξspec_constraints[1:np]

    ξ_constraints .= zbulk
    for j in 1:np
        ξj = viewn(ξ,nc,j)
        βj = β[j]
        ∑ξj = sum(ξj)
        β_constraints[j] = min(βj,1 - ∑ξj)
        for i in 1:nc
            ξ_constraints[i] -= βj*ξj[i]
        end
    end

    #fill specification constraints
    val1,spec1,val2,spec2 = state.val1,state.spec1,state.val2,state.spec2

    if spec1 == temperature
        output[end-1] = (T - val1)/val1
    elseif spec1 == pressure
        output[end-1] = (p_end - val1)/ps
    elseif spec1 == volume
        output[end-1] = (∑v - val1)/vs
    elseif spec1 isa Vfrac
        i_spec = spec1.k
        output[end-1] = β[i_spec] - val1
    else #general caloric term, valid for enthalpy , entropy, internal energy, gibbs, helmholtz
        output[end-1] = (val1 - ∑a - p_end*∑v - T*∑s)*RTinv
    end

    if spec2 == temperature
        output[end] = (T - val2)/val2
    elseif spec2 == pressure
        output[end] = (p_end - val2)/ps
    elseif spec2 == volume
        output[end] = (∑v - val2)/vs
    elseif spec2 isa Vfrac
        i_spec = spec2.k
        output[end-1] = β[i_spec] - val2
    else
        output[end] = (val2 - ∑a - p_end*∑v - T*∑s)*RTinv
    end
    #@show primalval(p_constraints)
    #@show primalval(μ_constraints)
    #@show primalval(β_constraints)
    #@show primalval(ξ_constraints)
    return output
end

function VT_chemical_potential_res_dual(model,V,T,z)
    TT = Base.promote_eltype(model,V,T,z)
    tag = ForwardDiff.Tag{ZVar{typeof(eos_res), typeof(model), TT, TT},eltype(z)}()
    N = Solvers.chunksize(z)
    return ForwardDiff.Dual{typeof(tag),TT,N}
end

function xy_flash_config(model::F1,input::F2) where {F1,F2}
    nc = length(model)
    chunk = ForwardDiff.Chunk(length(model))
    jconfig = ForwardDiff.JacobianConfig(XYFlashTag(),input,input,chunk)
    Solvers.chunksize(viewn(input,nc,1))
    ∂type = VT_chemical_potential_res_dual(model,input[1],input[1],viewn(jconfig.duals[1],nc,1))
    ∂x = similar(viewn(input,nc,1),∂type)
    gconfig = _gradient_config(∂x)
    return jconfig,gconfig
end

function _gradient_config(duals::AbstractVector{<:ForwardDiff.Dual{T,V,N}}) where {T,V,N}
    seeds = ForwardDiff.construct_seeds(ForwardDiff.Partials{N,V})
    ForwardDiff.GradientConfig{T,V,N,typeof(duals)}(seeds, duals)
end

function xy_flash(model::EoSModel,spec::FlashSpecifications,z,flash::FlashResult,method::FlashMethod)
    if numphases(flash) == 1
        throw(ArgumentError("xy_flash cannot use single phase initial points as starting points."))
    end
    comps0 = flash.compositions
    β0 = flash.fractions
    volumes0 = flash.volumes
    data0 = flash.data
    T0 = data0.T
    rtol = method.rtol
    atol = method.atol
    max_iters = method.max_iters
    return xy_flash(model,spec,z,comps0,β0,volumes0,T0;rtol,atol,max_iters)
end

function positive_linesearch(v, δ; τ = 1.0, decay = 0.5, tol = 1e-10)
    α = 1.0
    done = false
    while !done
        done = true
        for i in 1:length(v)
            vi,δi = v[i],δ[i]
            if vi + α * δi < (1 - τ) * vi
                if α < tol
                    return NaN
                end
                α *= decay
                done = false
            end
        end
    end
    return α
end

function minimize_df_linesearch(f!,F,v,vcache,Fnorm,α,decay = 0.5)
    done = false
    vcache .= α .* v
    while !done
        f!(F,vcache)
        Fnorm_i = norm(F,2)
        if Fnorm_i > Fnorm
            if α < 1e-4
                return NaN
            end
            α *= decay
            vcache .= α .* v
            done = false
        else
            done = true
        end
    end
    return α
end



function xy_flash(model::EoSModel,spec::FlashSpecifications,z,flash::FlashResult)
    if numphases(flash) == 1
        throw(ArgumentError("xy_flash cannot use single phase initial points as starting points."))
    end
    comps0 = flash.compositions
    β0 = flash.fractions
    volumes0 = flash.volumes
    T0 = flash.data.T
    return xy_flash(model,spec,z,comps0,β0,volumes0,T0)
end

function xy_flash(model::EoSModel,spec::FlashSpecifications,z,comps0,β0,volumes0,T0;rtol = 1e-12,atol = 1e-10,max_iters = 50)
    if length(comps0) == 1
        throw(ArgumentError("xy_flash cannot use single phase initial points as starting points."))
    end
    ∑z = sum(z)
    val1,val2 = spec.val1,spec.val2
    np = length(volumes0)
    nc = length(model)
    l = np*nc + 2*np + 1
    input = fill(zero(Base.promote_eltype(model,val1,val2,z)),l)

    #which inputs follow 0 < xi < 1
    input_0_1 = fill(false,l)
    idx_comps = 1:(nc*np)
    idx_volumes = 1:np .+ idx_comps[end]
    idx_β = 1:(np) .+ idx_volumes[end]
    input_0_1[idx_β] .= true
    input_0_1[idx_comps] .= true
    #we want to select the anchor phase with biggest fraction
    idx = sortperm(β0)
    snorm = Inf*one(eltype(input))
    flash_vars = xy_input_to_flash_vars(input,np,nc,z)
    compsx,βx,volumesx = flash_vars
    for j in 1:np
        k = idx[j]
        volumesx[j] = volumes0[k]
        βx[j] = β0[k]
        wxj = viewn(compsx,nc,j)
        wj0 = comps0[k]
        for i in 1:nc
            wxj[i] = wj0[i]
        end
        wxj ./= sum(wxj)
    end
    βx ./= sum(βx)
    if spec.spec1 == temperature 
        input[end] = spec.val1
    elseif spec.spec2 == temperature
        input[end] = spec.val2
    else
        input[end] = T0
    end
    _1 = one(eltype(input))
    new_spec = normalize_spec(set_vfrac(spec,idx),∑z*_1)
    zz = z * (1/∑z*_1)
    J = similar(input,(l,l))
    J .= 0
    piv = zeros(Int,l)
    x = input
    x_old = copy(input)
    x_cache = copy(input)
    F = copy(input)
    s = copy(input)
    #config,μconfig = xy_flash_config(model,input)
    f!(output,input) = xy_flash_neq(output,model,zz,np,input,new_spec,nothing)
    config = ForwardDiff.JacobianConfig(f!,F,x)
    #ForwardDiff.jacobian!(J,f!,F,x,config,Val{false}())
    f!(F,x)
    converged = false
    nan_converged = !all(isfinite,x)
    Fnorm = norm(F,Inf)
    Fnorm < rtol && (converged = true)

    for i in 1:max_iters
        converged && break
        nan_converged && break
        ForwardDiff.jacobian!(J,f!,F,x,config,Val{false}())
        lu = Solvers.unsafe_LU!(J,piv)
        s .= -F
        ldiv!(lu,s) #s .= J\F
        x_old .= x
        #bound 0-1 variables.
        
        α = positive_linesearch(@view(x[input_0_1]),@view(s[input_0_1]))
        #@show α
        #@show s
        #α = minimize_df_linesearch(f!,F,x,x_cache,norm(F,2),α)
        
        x .= x_old .+ α .* s
        snorm_old = snorm
        snorm = α*norm(s,2)
        if abs(snorm-snorm_old) < 1e-8
            converged = true
        end
        Fnorm = norm(F,Inf)
        
        xmax = maximum(x)
        xnorm = Solvers.dnorm(x,x_old,Inf)
        Fnorm < rtol && (converged = true)
        xnorm < atol && (converged = true)
        nan_converged = !all(isfinite,x)
    end

    comps_result,β_result,volumes_result,T_result = xy_input_to_result(x,np,nc,z)
    sp1,sp2 = spec.spec1,spec.spec2
    if sp1 == pressure
        p_result = spec.val1
    elseif sp2 == pressure
        p_result = spec.val2
    else
        p_result = pressure(model,volumes_result[end],T_result,comps_result[end])
    end
    β_result .*= ∑z
    return FlashResult(model,p_result,T_result,comps_result,β_result,volumes_result)
end

#======================

high level interface:

=======================#
"""
    GeneralizedXYFlash{T}(;kwargs...)

Method to solve non-reactive multicomponent, two-phase flash problem, using a generalized formulation.

Only two phases are supported. if `K0` is `nothing`, it will be calculated via fugacity coefficients at p,T conditions.

### Keyword Arguments:
- `equilibrium` (optional) = equilibrium type ":vle" for liquid vapor equilibria, ":lle" for liquid liquid equilibria, `:unknown` if not specified
- `p0` (optional), initial guess pressure, ignored if pressure is one of the flash specifications. 
- `T0` (optional), initial guess temperature, ignored if temperature is one of the flash specifications. 
- `K0` (optional), initial guess for the constants K
- `x0` (optional), initial guess for the composition of phase x
- `y0` = optional, initial guess for the composition of phase y
- `vol0` = optional, initial guesses for phase x and phase y volumes
- `atol` = absolute tolerance to stop the calculation
- `rtol` = relative tolerance to stop the calculation
- `max_iters` = maximum number of iterations
"""
struct GeneralizedXYFlash{P,T} <: FlashMethod
    equilibrium::Symbol
    T0::Union{P,Nothing}
    p0::Union{P,Nothing}
    K0::Union{Vector{T},Nothing}
    x0::Union{Vector{T},Nothing}
    y0::Union{Vector{T},Nothing}
    v0::Union{Tuple{T,T},Nothing}
    atol::Float64
    rtol::Float64
    max_iters::Int
end

Base.eltype(method::GeneralizedXYFlash{T}) where T = T

function index_reduction(m::GeneralizedXYFlash,idx::AbstractVector)
    equilibrium,T0,p0,K0,x0,y0,v0,atol,rtol,max_iters = m.equilibrium,m.T0,m.p0,m.K0,m.x0,m.y0,m.v0,m.atol,m.rtol,m.max_iters
    K0 !== nothing && (K0 = K0[idx])
    x0 !== nothing && (x0 = x0[idx])
    y0 !== nothing && (y0 = y0[idx])
    return GeneralizedXYFlash(;equilibrium,T0,p0,K0,x0,y0,v0,atol,rtol,max_iters)
end

numphases(::GeneralizedXYFlash) = 2

function GeneralizedXYFlash(;equilibrium = :unknown,
                        T0 = nothing,
                        p0 = nothing,
                        K0 = nothing,
                        x0 = nothing,
                        y0 = nothing,
                        v0 = nothing,
                        rtol = 1e-12,
                        atol = 1e-10,
                        max_iters = 100,
                        flash_result = nothing)
    !(is_vle(equilibrium) | is_lle(equilibrium) | is_unknown(equilibrium))  && throw(error("invalid equilibrium specification for GeneralizedXYFlash"))
    if flash_result isa FlashResult
        comps,β,volumes = flash_result.compositions,flash_result.fractions,flash_result.volumes
        @assert length(comps) == 2
        w1,w2 = comps[1],comps[2]
        v = (volumes[1],volumes[2])
        P00 = flash_result.data.p
        T00 = flash_result.data.T
        return GeneralizedXYFlash(;equilibrium = equilibrium,T0 = T00,p0 = P00,x0 = w1,y0 = w2,v0 = v,rtol = rtol,atol = atol,max_iters = max_iters)
    end
    
    if K0 == x0 == y0 === v0 == nothing #nothing specified
        #is_lle(equilibrium)
        T = Nothing
    else
        if !isnothing(K0) & isnothing(x0) & isnothing(y0) #K0 specified
            T = eltype(K0)
        elseif isnothing(K0) & !isnothing(x0) & !isnothing(y0)  #x0, y0 specified
            T = eltype(x0)
        else
            throw(error("invalid specification of initial points"))
        end
    end

    if T0 == nothing && p0 == nothing
        S = Nothing
    elseif T0 != nothing && p0 != nothing
        S = typeof(T0*p0)
    else
        S = typeof(something(T0,p0))
    end
    return GeneralizedXYFlash{S,T}(equilibrium,T0,p0,K0,x0,y0,v0,atol,rtol,max_iters)
end

function px_flash_x0(model,p,x,z,spec::F,method::GeneralizedXYFlash) where F
    if method.T0 == nothing
        T,_phase = _Tproperty(model,p,x,z,spec)
    else
        T = method.T0
        _phase = :eq #we suppose this
    end

    TT = Base.promote_eltype(model,p,x,z,T)
    if _phase != :eq
        return FlashResult(model,p,T,z,phase = _phase)
    end

    return pt_flash_x0(model,p,T,z,method;k0 = :suggest) 
end

function px_flash_pure(model,p,x,z,spec::F,T0 = nothing) where F
    Ts,vl,vv = saturation_temperature(model,p)
    ∑z = sum(z)
    x1 = SA[1.0]
    spec_to_vt(model,vl,T,x1,spec)
    xl = ∑z*spec_to_vt(model,vl,T,x1,spec)
    xv = ∑z*spec_to_vt(model,vv,T,x1,spec)
    βv = (x - xl)/(xv - xl)
    if !isfinite(βv)
        return FlashResultInvalid(1,βv)
    elseif βv < 0 || βv > 1
        phase0 = βv < 0 ? :liquid : :vapour
        T,_phase = _Tproperty(model,p,h/∑z,SA[1.0],spec,T0 = T0,phase = phase0)
        return FlashResult(model,p,T,∑z,phase = _phase)
    else
        build_flash_result_pure(model,p,Ts,z,vl,vv,βv)
    end
end

function tx_flash_x0(model,T,x,z,spec::F,method::GeneralizedXYFlash) where F
    if method.p0 == nothing
        p,_phase = _Pproperty(model,T,x,z,spec)
    else
        p = method.p0
        _phase = :eq #we suppose this
    end
    
    TT = Base.promote_eltype(model,T,x,z,T)
    if _phase != :eq
        return FlashResult(model,p,T,z,phase = _phase)
    end

    return pt_flash_x0(model,p,T,z,method;k0 = :suggest)
end

function tx_flash_pure(model,T,x,z,spec::F,P0 = nothing) where F
    ps,vl,vv = saturation_pressure(model,T)
    ∑z = sum(z)
    x1 = SA[1.0]
    spec_to_vt(model,vl,T,x1,spec)
    xl = ∑z*spec_to_vt(model,vl,T,x1,spec)
    xv = ∑z*spec_to_vt(model,vv,T,x1,spec)
    βv = (x - xl)/(xv - xl)
    if !isfinite(βv)
        return FlashResultInvalid(1,βv)
    elseif βv < 0 || βv > 1
        phase0 = βv < 0 ? :liquid : :vapour
        p,_phase = _Pproperty(model,T,x/∑z,SA[1.0],spec,p0 = P0,phase = phase0)
        return FlashResult(model,p,T,∑z,phase = _phase)
    else
        build_flash_result_pure(model,p,Ts,z,vl,vv,βv)
    end
end

function βflash_pure(model,spec::F,x,βv,z) where F
    if spec == pressure
        T,vl,vv = saturation_temperature(model,x)
        p = x
    elseif spec == temperature
        p,vl,vv = saturation_pressure(model,T)
        T = x
    else
        throw(ArgumentError("invalid specification for βflash_pure: $spec"))
    end
    ∑z = sum(z)
    x1 = SA[1.0]
    if !isfinite(βv)
        return FlashResultInvalid(1,βv)
    elseif βv == 1
        return FlashResult([[1.0]],[∑z],[vv],FlashData(p,T))
    elseif βv == 0
        return FlashResult([[1.0]],[∑z],[vl],FlashData(p,T))   
    elseif βv < 0 || βv > 1
        throw(error("invalid specification of vapour fraction, it must be between 0 and 1."))
    else
        build_flash_result_pure(model,p,T,z,vl,vv,∑z*βv)
    end
end