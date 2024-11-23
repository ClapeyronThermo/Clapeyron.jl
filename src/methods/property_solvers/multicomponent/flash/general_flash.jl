struct FlashSpecifications{S1,V1,S2,V2}
    spec1::S1
    val1::V1
    spec2::S2
    val2::V2

    function FlashSpecifications(spec1::S1,val1::V1,spec2::S2,val2::V2) where {S1,V1,S2,V2}
        @assert spec1 == volume || spec1 == temperature || spec1 == entropy || spec1 == enthalpy || spec1 == internal_energy || spec1 == pressure
        @assert spec2 == volume || spec2 == temperature || spec2 == entropy || spec2 == enthalpy || spec2 == internal_energy || spec2 == pressure
        @assert spec1 !== spec2
        return new{S1,V1,S2,V2}(spec1,val1,spec2,val2)
    end
end

function normalize_spec(spec::FlashSpecifications,k)
    return FlashSpecifications(spec.spec1,spec.val1*k,spec.spec2,spec.val2*k)
end

function xy_input_to_flash_vars(input,np,nc,z)
    idx_comps_end = (nc - 1)*(np - 1)
    idx_comps = 1:idx_comps_end
    idx_volumes = (1:np) .+ idx_comps_end
    idx_β = (1:(np-1)) .+ idx_volumes[end]
    comps0 = @view input[idx_comps]
    volumes = @view input[idx_volumes]
    β = @view input[idx_β]
    #fill last component vector
    comps0_end = similar(input,nc-1)
    comps0_end .= @view z[1:end-1]
    β_end = 1 - sum(β)
    for i in 1:(nc - 1)
        _xi = z[i]
        
        res = zero(eltype(input))
        res += z[i]
        for j in 1:(np - 1)
            compsj = viewn(comps0,nc-1,j)
            res -= β[j]*compsj[i]
        end
        comps0_end[i] = res
    end

    #=
    for j in 1:(np - 1)
        compsj = viewn(comps0,nc-1,j)
        ∑xj_end = 1 - sum(compsj)
        for i in 1:(nc - 1)
            comps0_end[i] -= β_end*∑xj_end
        end
    end =#

    #clamp negative values
    for i in 1:nc-1
        if comps0_end[i] < 0
            comps0_end[i] = 0
        end
    end

    return comps0,β,volumes,comps0_end,β_end
end

function xy_input_to_result(input,np,nc,z)
    compsx,βx,volumes,_xi,β_end = xy_input_to_flash_vars(input,np,nc,z)
    TT = eltype(input)

    comps = Vector{Vector{TT}}(undef,np)
    β = Vector{TT}(undef,np)
    for j in 1:np-1
        compsi = viewn(compsx,nc-1,j)
        comps[j] = collect(FractionVector(compsi))
        β[j] = βx[j]
    end
    comps[end] = collect(FractionVector(_xi))
    β[end] = β_end
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

function xy_flash_neq(output,model,zbulk,np,input,state::F) where F
    #=
    variables:

    phase fractions (β) [np]
    phase compositions (xi) [np*nc]
    phase volumes (vi) [np]
    T [1]
    total: np*nc + 2*np + 1
    ----------------------------------------
    equations:
    pressure equalities: [np - 1]
    chemical potential equalities: [(np-1)*nc]
    enthalpy equal to especification: [1]
    pressure equal to specification: [1]

    total: np*nc + np - nc + 1
    ----------------------------------------
    missing: np + nc
    ----------------------------------------
    simplifications:
    x_end = (z - sum(βi*xi))/β_end [-nc]
    β_end = 1 - sum(βi) [-1]
    sum(xi) = 1 [-(np - 1)]


    restrictions (to add in solver):
    0 < βi < 1
    0 <. xi .< 1
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

    =#
    T = input[end]
    nc = length(model)

    flash_vars = xy_input_to_flash_vars(input,np,nc,zbulk)
    comps0,β,volumes,comps0_end,β_end = flash_vars
    w_end = FractionVector(comps0_end)
    v_end = volumes[end]
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

    ∑a = β_end*a_end
    if needs_pv
        ∑v = β_end*v_end
    else
        ∑v = zero(∑a)
    end
    ∑s = β_end*s_end
    ps = p_scale(model,zbulk)
    R = Rgas(model)
    RT = R*T
    RTinv = 1/RT
    vs = RT/ps
    #fill pressure constraints:
    idx_p_constraints = 1:(np-1)
    p_constraints = @view output[idx_p_constraints]
    for j in 1:(np-1)
        βj = β[j]
        vj = volumes[j]
        if needs_pv
            ∑v += βj*vj
        end
        comps0j = viewn(comps0,nc-1,j)
        wj = FractionVector(comps0j)
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
        ∑a += βj*aj
        ∑s += βj*sj
        p_constraints[j] = (pj - p_end)/ps
    end

    #fill chemical potential equalities:

    idx_μ_constraints = (1:((np-1)*nc)) .+ (np - 1)
    μ_constraints = @view output[idx_μ_constraints]
    μ_end = similar(output,nc)

    VT_chemical_potential_res!(μ_end,model,v_end,T,FractionVector(comps0_end))

    for j in 1:(np - 1)
        Fj = @inbounds viewn(μ_constraints,nc,j)
        for i in 1:nc
            Fj[i] = μ_end[i]
        end
    end

    for j in 1:(np - 1)
        vj = @inbounds volumes[j]
        comps0j = viewn(comps0,nc-1,j)
        wj = FractionVector(comps0j)
        VT_chemical_potential_res!(μ_end,model,vj,T,wj)
        Fj = viewn(μ_constraints,nc,j)
        for i in 1:nc
            μ1i = Fj[i]
            μji = μ_end[i]
            Δuᵣ = μ1i - μji
            Δu = Δuᵣ*RTinv + log(vj*w_end[i]/(v_end*wj[i]))
            Fj[i] = Δu
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
    else #general caloric term, valid for enthalpy , entropy, internal energy, gibbs, helmholtz
        output[end-1] = (val1 - ∑a - p_end*∑v - T*∑s)*RTinv
    end

    if spec2 == temperature
        output[end] = (T - val2)/val2
    elseif spec2 == pressure
        output[end] = (p - val2)/ps
    elseif spec2 == volume
        output[end] = (∑v - val2)/vs
    else
        output[end] = (val2 - ∑a - p_end*∑v - T*∑s)*RTinv
    end

    return output
end

function xy_flash(model::EoSModel,spec::FlashSpecifications,z,flash::FlashResult,method::FlashMethod)
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

function xy_flash(model::EoSModel,spec::FlashSpecifications,z,comps0,β0,volumes0,T0;rtol = 1e-12,atol = 1e-10,max_iters = 50)
    ∑z = sum(z)
    val1,val2 = spec.val1,spec.val2
    np = length(volumes0)
    nc = length(model)
    l = np*nc + np - nc + 1
    input = fill(zero(Base.promote_eltype(model,val1,val2,z)),l)

    #which inputs follow 0 < xi < 1
    input_0_1 = fill(false,l)
    idx_comps_end = (nc - 1)*(np - 1)
    idx_comps = 1:idx_comps_end
    idx_volumes = 1:np .+ idx_comps_end
    idx_β = 1:(np-1) .+ idx_volumes[end]
    input_0_1[idx_β] .= true
    input_0_1[idx_comps] .= true
    #we want to select the anchor phase with biggest fraction
    idx = sortperm(β0)
    flash_vars = xy_input_to_flash_vars(input,np,nc,z)
    compsx,βx,volumesx,_,_ = flash_vars
    for j in 1:(np-1)
        k = idx[j]
        volumesx[j] = volumes0[k]
        βx[j] = β0[k]/∑z
        wxj = viewn(compsx,nc-1,j)
        wj0 = comps0[k]
        for i in 1:nc-1
            wxj[i] = wj0[i]
        end
    end
    volumesx[end] = volumes0[idx[end]]

    if spec.spec1 == temperature 
        input[end] = spec.val1
    elseif spec.spec2 == temperature
        input[end] = spec.val2
    else
        input[end] = T0
    end
    
    _1 = one(eltype(input))
    is_caloric = requires_st(spec) || requires_a(spec)
    if is_caloric
        new_spec = normalize_spec(spec,1/∑z*_1)
    else
        new_spec = normalize_spec(spec,_1)
    end

    f!(output,input) = xy_flash_neq(output,model,z,np,input,new_spec)
    J = similar(input,(l,l))
    J .= 0
    piv = zeros(Int,l)
    x = input
    x_old = copy(input)
    F = copy(input)
    s = copy(input)
    f!(F,x)
    converged = false
    nan_converged = !all(isfinite,x)
    if norm(F,Inf) < rtol
        converged = true
    end
    config = ForwardDiff.JacobianConfig(f!,F,x)
    for i in 1:max_iters
        converged && break
        nan_converged && break
        ForwardDiff.jacobian!(J,f!,F,x,config)
        lu = Solvers.unsafe_LU!(J,piv)
        s .= F
        ldiv!(lu,s) #s .= J\F
        x_old .= x
        #TODO: backtracking? trust region?, reformulate in SS form if the step is too far?  
        x .= x_old .- s
        #bound 0-1 variables.
        for j in 1:l
            check_bounds = input_0_1[j]
            if check_bounds && x[j] < 0
                x[j] = 0.5*x_old[j]
            elseif check_bounds && x[j] < 0
                x[j] = 0.5*(1 + x_old[j])
            end
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
        p_result = pressure(model,volumes_result[end],T,comps_result[end])
    end
    β_result .*= ∑z
    return FlashResult(model,p_result,T_result,comps_result,β_result,volumes_result)
end

#==

high level interface

==#

"""
    GeneralizedPXFlash{T}(;kwargs...)

Method to solve non-reactive multicomponent flash problem, using a generalized formulation.

Only two phases are supported. if `K0` is `nothing`, it will be calculated via fugacity coefficients at p,T conditions.

### Keyword Arguments:
- equilibrium = equilibrium type ":vle" for liquid vapor equilibria, ":lle" for liquid liquid equilibria
- `T0` (optional), initial guess temperature
- `K0` (optional), initial guess for the constants K
- `x0` (optional), initial guess for the composition of phase x
- `y0` = optional, initial guess for the composition of phase y
- `vol0` = optional, initial guesses for phase x and phase y volumes
- `atol` = absolute tolerance to stop the calculation
- `rtol` = relative tolerance to stop the calculation
- `max_iters` = maximum number of iterations
"""
struct GeneralizedPXFlash{S,T} <: FlashMethod
    equilibrium::Symbol
    T0::Union{S,Nothing}
    K0::Union{Vector{T},Nothing}
    x0::Union{Vector{T},Nothing}
    y0::Union{Vector{T},Nothing}
    v0::Union{Tuple{T,T},Nothing}
    atol::Float64
    rtol::Float64
    max_iters::Int
end

Base.eltype(method::GeneralizedPXFlash{T}) where T = T

function index_reduction(m::GeneralizedPXFlash,idx::AbstractVector)
    equilibrium,T0,K0,x0,y0,v0,atol,rtol,max_iters = m.equilibrium,m.T0,m.K0,m.x0,m.y0,m.v0,m.atol,m.rtol,m.max_iters
    K0 !== nothing && (K0 = K0[idx])
    x0 !== nothing && (x0 = x0[idx])
    y0 !== nothing && (y0 = y0[idx])
    return GeneralizedPXFlash(;equilibrium,T0,K0,x0,y0,v0,atol,rtol,max_iters)
end

numphases(::GeneralizedPXFlash) = 2

function GeneralizedPXFlash(;equilibrium = :unknown,
                        T0 = nothing,
                        K0 = nothing,
                        x0 = nothing,
                        y0 = nothing,
                        v0 = nothing,
                        rtol = 1e-12,
                        atol = 1e-10,
                        max_iters = 100)
    !(is_vle(equilibrium) | is_lle(equilibrium) | is_unknown(equilibrium))  && throw(error("invalid equilibrium specification for GeneralizedPXFlash"))
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
    S = typeof(T0)
    return GeneralizedPXFlash{S,T}(equilibrium,T0,K0,x0,y0,v0,atol,rtol,max_iters)
end

is_vle(method::GeneralizedPXFlash) = is_vle(method.equilibrium)
is_lle(method::GeneralizedPXFlash) = is_lle(method.equilibrium)

function px_flash_x0(model,p,x,z,spec::F,method::GeneralizedPXFlash) where F
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
    if (0 <= βv <= 1)
    elseif βv > 1
        return build_flash_result_pure(model,p,Ts,z,vl,vv,βv)
    elseif !isfinite(βv)
        return FlashResultInvalid(1,βv)
    elseif βv < 0 || βv > 1
        phase0 = βv < 0 ? :liquid : :vapour
        T,_phase = _Tproperty(model,p,h/∑z,SA[1.0],spec,T0 = T0,phase = phase0)
        return FlashResult(model,p,T,∑z,phase = _phase)
    else
        build_flash_result_pure(model,p,Ts,z,vl,vv,βv)
    end
end