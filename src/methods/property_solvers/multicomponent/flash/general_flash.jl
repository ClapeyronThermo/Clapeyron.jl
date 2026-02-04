struct Vfrac
    k::Int
end

"""
    FlashSpecifications(;v,T,p,h,s,u,q)

Struct that holds two specifications for a general flash.
The keyword arguments have the following meaning:

- `T`: temperature `[K]`
- `v`: total volume `[m³]`
- `p`: pressure `[Pa]`
- `h`: enthalpy `[J]`
- `s`: entropy `[J K⁻¹]`
- `u`: internal energy `[J]`
- `q`: vapour fraction

## Examples:
```
specs = FlashSpecifications(p = 101325,T = 298.15) #PT flash
specs = FlashSpecifications(Clapeyron.pressure,101325,Clapeyron.temperature,298.15) #equivalent
```
"""
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
spec_intensive(::Vfrac) = true

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

function xy_input_to_flash_vars(input,np,nc,comps_offset = 0)
    idx_comps_end = np*(nc-comps_offset)

    idx_comps = 1:idx_comps_end
    idx_volumes = (1:np) .+ idx_comps[end]
    idx_β = (1:np) .+ idx_volumes[end]

    comps = @view input[idx_comps]
    volumes = @view input[idx_volumes]
    β = @view input[idx_β]

    return comps,β,volumes
end

function xy_input_to_result(spec,input,np,nc,z)
    compsx,βx,volumes = xy_input_to_flash_vars(input,np,nc)
    TT = eltype(input)
    comps = Vector{Vector{TT}}(undef,np)
    β = Vector{TT}(undef,np)
    for j in 1:np
        compsi = viewn(compsx,nc,j)
        ∑ξ = sum(compsi)
        comps[j] = collect(compsi)
        comps[j] ./= ∑ξ
        volumes[j] = volumes[j] #we already normalize inside xy_flash
        β[j] = βx[j]
    end
    Tres = input[end]

    if spec.spec1 == temperature
        T = spec.val1
    elseif spec.spec2 == temperature
        T = spec.val2
    else
        T = Tres
    end
    #we return in order of increasing molar volumes
    idx = sortperm(volumes)
    return comps[idx],β[idx],volumes[idx],T
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

function __min(x,y)
    if x < y
        return x
    elseif y < x
        return y
    else
        return y
    end
end

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

    flash_vars = xy_input_to_flash_vars(input,np,nc)
    ξ,β,volumes = flash_vars
    #create fractions
    frac = similar(ξ)
    for j in 1:np
        xj,ξj = viewn(frac,nc,j),viewn(ξ,nc,j)
        xj .= ξj
        xj ./= sum(xj)
    end
    #@show primalval(ξ)
    #@show primalval(β)
    #@show primalval(volumes)
    β_end,j_end = β[np],np
    w_end = viewn(ξ,nc,j_end)
    x_end = viewn(frac,nc,j_end)
    v_end = volumes[j_end]
    needs_a = requires_a(state)
    needs_st = requires_st(state)
    needs_pv = requires_pv(state)

    if needs_st
        a_end,dadv_end,dadT_end = ∂f_vec(model,v_end,T,x_end)
        p_end,s_end = -dadv_end,-dadT_end
    elseif needs_a && !needs_st
        a_end,dadv_end = f∂fdV(model,v_end,T,x_end)
        p_end,s_end = -dadv_end,zero(a_end)
    else
        p_end = pressure(model,v_end,T,x_end)
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
        wj = viewn(ξ,nc,jj)
        xj = viewn(frac,nc,jj)
        vj = volumes[j]
        if needs_pv
            ∑v += βj*vj
        end

        if needs_st
            aj,dadvj,dadTj = ∂f_vec(model,vj,T,xj)
            pj,sj = -dadvj,-dadTj
        elseif needs_a && !needs_st
            aj,dadvj = f∂fdV(model,vj,T,xj)
            pj,sj = -dadvj,zero(aj)
        else
            pj = pressure(model,vj,T,xj)
            sj,aj = zero(pj),zero(pj)
        end
        #@show primalval(pj)
        ∑a += βj*aj
        ∑s += βj*sj
        p_constraints[jj] = (pj - p_end)/ps
    end
    #fill chemical potential equalities:

    idx_μ_constraints = (1:((np-1)*nc)) .+ (idx_p_constraints[end])
    μ_constraints = @view output[idx_μ_constraints]
    μ_end = similar(output,nc)

    VT_chemical_potential_res!(μ_end,model,v_end,T,x_end)

    jj = 0
    for j in 1:np
        j == j_end && continue
        jj += 1
        Fj = @inbounds viewn(μ_constraints,nc,jj)
        for i in 1:nc
            #Fj[i] = exp(μ_end[i]*RTinv)*w_end[i]*vs/v_end
            Fj[i] = μ_end[i]
        end
    end
    jj = 0
    for j in 1:np
        j == j_end && continue
        jj += 1
        vj = @inbounds volumes[j]
        wj = viewn(ξ,nc,j)
        xj = viewn(frac,nc,j)
        VT_chemical_potential_res!(μ_end,model,vj,T,xj)
        Fj = viewn(μ_constraints,nc,jj)
        for i in 1:nc
            μ1i = Fj[i]
            μji = μ_end[i]
            Δuᵣ = μ1i - μji
            Δu = Δuᵣ*RTinv + log(vj) + log(x_end[i]) - log(v_end) - log(xj[i])
            Fj[i] = Δu
        end
    end

    #fill β,extended composition constraints
    idx_β_constraints = (1:np) .+ (idx_μ_constraints[end])
    idx_ξ_constraints = (1:nc) .+ (idx_β_constraints[end])
    
    β_constraints = @view output[idx_β_constraints]
    ξ_constraints = @view output[idx_ξ_constraints]

    ξ_constraints .= zbulk

    val1,spec1,val2,spec2 = state.val1,state.spec1,state.val2,state.spec2

    if spec1 isa Vfrac
        idx_β = spec1.k
        βx = val1
    elseif spec2 isa Vfrac
        idx_β = spec2.k
    else
        idx_β = 0
        βx = val2
    end

    for j in 1:np
        ξj = viewn(ξ,nc,j)
        βj = β[j]
        ∑ξj = sum(ξj)
        #TODO: using min(a,b) in conjunction with AD seems to be equivalent in the Newton-min method, a way to solve MCP.
        #there are better ways to solve MCP.
        #ben-gharbia uses a Non-Parametric-Interior-Point method (npipm)
        #there is also MixedComplementarityProblems.jl
        #@show primalval(βj),primalval(1 - ∑ξj)
        β_constraints[j] = __min(βj,1 - ∑ξj)

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
        #output[end-1] = (∑v - val1)/vs
        output[end-1] = log(∑v/val1)  #this was the best residual for volume
        #output[end-1] = vs/∑v - vs/val1
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
        #output[end] = (∑v - val2)/vs
        output[end] = log(∑v/val2) #this was the best residual for volume
        #output[end] = vs/∑v - vs/val2
    elseif spec2 isa Vfrac
        i_spec = spec2.k
        output[end-1] = β[i_spec] - val2
    else
        output[end] = (val2 - ∑a - p_end*∑v - T*∑s)*RTinv
    end
    return output
end

#the idea is to not update the x at specification values
function detect_and_set_slack_variables!(x,spec::FlashSpecifications,np,nc,comps_offset = 0)
    slack = similar(x,Bool)
    slack .= false
    if spec.spec1 == temperature
        slack[end] = true
        x[end] = spec.val1
    elseif spec.spec2 == temperature
        slack[end] = true
        x[end] = spec.val2
    end
    return slack
    #FIXME: we need to perform dew temperatures correctly, and that will need extra slacks.
    slack_comps = @view slack[1:np*(nc-comps_offset)]

    #bubbledew condition in first spec
    if spec.spec1 isa Vfrac && np == 2 && (iszero(spec.val1) || isone(spec.val1))
        k = spec.spec1.k
        xk = viewn(slack_comps,nc,k)
        xk .= true
    end

    #bubbledew condition in second spec
    if spec.spec2 isa Vfrac && np == 2 && (iszero(spec.val2) || isone(spec.val2))
        k = spec.spec2.k
        xk = viewn(slack_comps,nc,k)
        xk .= true
    end
    return slack
end

#we set the value of F[slack] = 0, J[slack_i,slack_i] = 1,  J[:,slack_i] = 0, J[slack_i,:] = 0


"""
    xy_flash(model,spec::FlashSpecifications,z,w0::FlashResult,method::GeneralizedXYFlash)
    xy_flash(model,spec::FlashSpecifications,z,w0::FlashResult;rtol = 1e-12,atol = 1e-10,max_iters = 50)

Routine to solve a non-reactive multiphase, multicomponent flash problem with arbitrary specifications.
Based on the unified formulation proposed by Ben Gharbia (2021).
The two necessary specifications are taken from `specs`, the initial point is obtained from `w0`.
Returns a [`FlashResult`](@ref) struct containing molar fractions, vapour fractions, molar volumes and the equilibrium temperature and pressure.

## Examples
```
spec = FlashSpecifications(p = 101325.0, T = 200.15) #p-T flash
model = cPR(["ethane","propane"],idealmodel=ReidIdeal)
z = [0.5,0.5] #bulk composition
x1 = [0.25,0.75] #liquid composition
x2 = [0.75,0.25] #gas composition
compositions = [x1,x2]
volumes = [6.44e-5,0.016]
fractions = [0.5,0.5]
p0,T0 = NaN,NaN #in p-T flash, pressure and temperature are already specifications
data = FlashData(p0,T0)
result0 = FlashResult(compositions,fractions,volumes,data) #a FlashResult containing all necessary information
result = xy_flash(model,spec,z,result0) #perform the flash
```
"""
function xy_flash end

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

function xy_flash(model::EoSModel,spec::FlashSpecifications,z,flash::FlashResult;rtol = 1e-14,atol = 1e-12,max_iters = 50)
    if numphases(flash) == 1
        throw(ArgumentError("xy_flash cannot use single phase initial points as starting points."))
    end
    comps0 = flash.compositions
    β0 = flash.fractions
    volumes0 = flash.volumes
    T0 = flash.data.T
    return xy_flash(model,spec,z,comps0,β0,volumes0,T0;rtol,atol,max_iters)
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
    flash_vars = xy_input_to_flash_vars(input,np,nc)
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
        wxj ./= sum(wj0)
    end
    βx ./= sum(βx)
    _1 = one(eltype(input))
    #normalize to 1 mol base
    new_spec = normalize_spec(set_vfrac(spec,idx),∑z*_1)
    #variables already set by the specifications
    slacks = detect_and_set_slack_variables!(input,spec,np,nc)
    !slacks[end] && (input[end] = T0)
    zz = z * (1/∑z*_1)

    J = similar(input,(l,l))
    J .= 0
    Jcache = similar(J)
    piv = zeros(Int,l)
    x = input
    x_old = copy(input)
    F = copy(input)
    s = copy(input)
    #config,μconfig = xy_flash_config(model,input)
    f!(output,input) = xy_flash_neq(output,model,zz,np,input,new_spec,nothing)
    srtol = abs2(cbrt(rtol))
    function Θ(_f,_z)
        f!(_f,_z)
        _f[slacks] .= 0
        0.5*dot(_f,_f)
    end
    config = ForwardDiff.JacobianConfig(f!,F,x)
    #ForwardDiff.jacobian!(J,f!,F,x,config,Val{false}())
    Θx = Θ(F,x)
    Fnorm = sqrt(2*Θx)
    spec_norm = norm(viewlast(F,2),Inf)
    converged = Fnorm < rtol
    nan_converged = !all(isfinite,x)
    max_iters_reached = false
    ii = 0
    for i in 1:max_iters
        ii = i
        converged && break
        nan_converged && break
        max_iters_reached && break
        ForwardDiff.jacobian!(J,f!,F,x,config,Val{false}())
        #TODO: fix volumes when they enter an unstable state. do it right here, were we have jacobian info.
        #do not iterate on slack variables
        Solvers.remove_slacks!(F,J,slacks)
        Jcache .= J
        #try to do LU, if it does not work, use modified SVD
        finite_F = all(isfinite,F)
        finite_J = all(isfinite,J)
        lu = Solvers.unsafe_LU!(J,piv)
        s .= -F
        ldiv!(lu,s)
        finite_s = all(isfinite,s)
        if !finite_s && finite_J
            JJ = svd(Jcache)
            S = JJ.S
            for i in eachindex(S)
                if S[i] == 0
                    S[i] = 1
                end
            end
            s .= -F
            ldiv!(JJ,s)
            s .= JJ\-F
        end

        for i in 1:length(s)
            si = s[i]
            abs(si) < eps(eltype(s)) && si < 0 &&  (s[i] = 0)
        end

        x_old .= x
        #bound positivity
        α0 = Solvers.positive_linesearch(x,s,decay = 0.8)
        #backtrack linesearch, so the next result is strictly better than the last
        α,Θx = Solvers.backtracking_linesearch!(Θ,F,x_old,s,Θx,x,α0,ignore = slacks)
        Fnorm = sqrt(2*Θx)
        spec_norm = norm(viewlast(F,2),Inf)
        snorm_old = snorm
        snorm = α*norm(s,2)
        δs = max(abs(snorm-snorm_old),spec_norm)
        Fnorm = norm(F,Inf)
        xnorm = Solvers.dnorm(x,x_old,Inf)
        converged = (Fnorm < rtol || xnorm < atol || δs < rtol)
        nan_converged = !all(isfinite,x) || !all(isfinite,s)
        isnan(δs) && (nan_converged = true)
        i == max_iters && (max_iters_reached = true)
    end
    if !converged || nan_converged || max_iters_reached
        x .= NaN
    end

    comps_result,β_result,volumes_result,T_result = xy_input_to_result(spec,x,np,nc,z)
    sp1,sp2 = spec.spec1,spec.spec2
    if sp1 == pressure
        p_result = spec.val1
    elseif sp2 == pressure
        p_result = spec.val2
    else
        p_result = pressure(model,volumes_result[end],T_result,comps_result[end])
    end
    β_result .*= ∑z
    flash_result = FlashResult(model,p_result,T_result,comps_result,β_result,volumes_result)
    return merge_duplicate_phases!(flash_result)
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
- `K0` (optional), initial guess for the K-values.
- `x0` (optional), initial guess for the composition of phase x.
- `y0` = optional, initial guess for the composition of phase y.
- `vol0` = optional, initial guesses for phase x and phase y volumes.
- `atol` = absolute tolerance to stop the calculation.
- `rtol` = relative tolerance to stop the calculation.
- `max_iters` = maximum number of iterations
- `flash_result::FlashResult`: can be provided instead of `x0`,`y0` and `vol0` for initial guesses.
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

index_reduction(m::GeneralizedXYFlash{Nothing,Nothing},idx::AbstractVector) = m

numphases(::GeneralizedXYFlash) = 2

function GeneralizedXYFlash(;equilibrium = :unknown,
                        T0 = nothing,
                        p0 = nothing,
                        K0 = nothing,
                        x0 = nothing,
                        y0 = nothing,
                        v0 = nothing,
                        rtol = 1e-14,
                        atol = 1e-12,
                        max_iters = 100,
                        flash_result = nothing)
    !(is_vle(equilibrium) | is_lle(equilibrium) | is_unknown(equilibrium))  && throw(error("invalid equilibrium specification for GeneralizedXYFlash"))
    if flash_result isa FlashResult
        comps,β,volumes = flash_result.compositions,flash_result.fractions,flash_result.volumes
        np = numphases(flash_result)
        np != 2 && incorrect_np_flash_error(GeneralizedXYFlash,flash_result)
        w1,w2 = comps[1],comps[2]
        v = (volumes[1],volumes[2])
        P00 = flash_result.data.p
        T00 = flash_result.data.T
        return GeneralizedXYFlash(;equilibrium = equilibrium,T0 = T00,p0 = P00,x0 = w1,y0 = w2,v0 = v,rtol = rtol,atol = atol,max_iters = max_iters)
    end

    if K0 == x0 == y0 === nothing #nothing specified
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

    if T == Nothing && v0 !== nothing
        TT = Base.promote_eltype(v0[1],v0[2])
        _v0 = (v0[1],v0[2])
    elseif T != nothing && v0 !== nothing
        TT = Base.promote_eltype(one(T),v0[1],v0[2])
        _v0 = (v0[1],v0[2])
    else
        TT = T
        _v0 = v0
    end

    if T0 === nothing && p0 === nothing
        S = Nothing
    elseif T0 !== nothing && p0 !== nothing
        S = typeof(T0*p0)
    else
        S = typeof(something(T0,p0))
    end
    return GeneralizedXYFlash{S,TT}(equilibrium,T0,p0,K0,x0,y0,_v0,atol,rtol,max_iters)
end

function px_flash_x0(model,p,x,z,spec::F,method::GeneralizedXYFlash) where F

    if spec == temperature
        T,_phase = x,:eq #we suppose equilibria
    elseif method.T0 === nothing
        T,_phase = _Tproperty(model,p,x,z,spec)
    else
        T,_phase = method.T0,:eq #we suppose equilibria
    end

    TT = Base.promote_eltype(model,p,x,z,T)
    if _phase != :eq
        return FlashResult(model,p,T,z,phase = _phase)
    end

    return pt_flash_x0(model,p,T,z,method)
end

function px_flash_pure(model,p,x,z,spec::F,T0 = nothing) where F

    ∑z = sum(z)
    x1 = SVector(1.0*one(∑z))
    TT = Base.promote_eltype(model,p,x,z)

    sat,crit,status = _extended_saturation_temperature(model,p)

    if status == :fail    
        return FlashResultInvalid(x1,one(TT))
    end

    if status == :supercritical
        Tc,Pc,Vc = crit
        if T0 !== nothing
            Tcrit0 = TT(T0)
        else
            Tcrit0 = TT(1.001Tc) #some eos have problems at exactly the critical point (SingleFluid("R123"))
        end
        Tsc,_phase = __Tproperty(model,p,x/∑z,x1,spec,:unknown,Tcrit0)
        return FlashResult(model,p,Tsc,SA[∑z*one(p)*one(Tsc)],phase = _phase)
    end

    Ts,vl,vv = sat

    xl = ∑z*spec_to_vt(model,vl,Ts,x1,spec)
    xv = ∑z*spec_to_vt(model,vv,Ts,x1,spec)
    βv = (x - xl)/(xv - xl)

    if !isfinite(βv)
        return FlashResultInvalid(x1,βv)
    elseif βv < 0 || βv > 1
        phase0 = βv < 0 ? :liquid : :vapour
        _T0 = T0 === nothing ? TT(Ts) : TT(primalval(T0))
        Tx,_phase = __Tproperty(model,p,x/∑z,x1,spec,phase0,_T0)
        return FlashResult(model,p,Tx,SA[∑z*one(p)*one(Tx)],phase = _phase)
    else
        return FlashResult(model,p,Ts,[x1,x1],[∑z-∑z*βv,∑z*βv],[vl,vv];sort = false)
    end
end

function tx_flash_x0(model,T,x,z,spec::F,method::GeneralizedXYFlash) where F

    if spec == pressure
        p,_phase = x,:eq #we suppose equilibria
    elseif method.p0 === nothing
        p,_phase = _Pproperty(model,T,x,z,spec)
    else
        p,_phase = x,:eq #we suppose equilibria
    end

    TT = Base.promote_eltype(model,T,x,z,T)
    if _phase != :eq
        return FlashResult(model,p,T,z,phase = _phase)
    end

    return pt_flash_x0(model,p,T,z,method)
end

function tx_flash_pure(model,T,x,z,spec::F,P0 = nothing) where F

    ∑z = sum(z)
    x1 = SA[1.0*one(∑z)]
    TT = Base.promote_eltype(model,T,x,z)

    sat,crit,status = _extended_saturation_pressure(model,T)

    if status == :fail
        return FlashResultInvalid(x1,one(TT))
    end

    if status == :supercritical
        Tc,Pc,Vc = crit #TODO: maybe use critical extrapolation instead?
        if P0 !== nothing
            Pcrit0 = TT(P0)
        else
            Pcrit0 = TT(1.001Pc) #some eos have problems at exactly the critical point (SingleFluid("R123"))
        end
        psc,_phase = __Pproperty(model,T,x/∑z,x1,spec,:unknown,Pcrit0)
        return FlashResult(model,psc,T,SA[∑z*one(psc)*one(T)])
    end

    ps,vl,vv = TT.(sat)

    xl = ∑z*spec_to_vt(model,vl,T,x1,spec)
    xv = ∑z*spec_to_vt(model,vv,T,x1,spec)
    βv = (x - xl)/(xv - xl)

    if !isfinite(βv)
        return FlashResultInvalid(x1,βv)
    elseif βv < 0 || βv > 1
        phase0 = βv < 0 ? :liquid : :vapour
        _p0 = P0 === nothing ? TT(ps) : TT(primalval(P0))
        px,_phase = __Pproperty(model,T,x/∑z,x1,spec,phase0,_p0)
        return FlashResult(model,px,T,SA[∑z*one(px)*one(T)],phase = _phase)
    else
        return FlashResult(model,ps,T,[x1,x1],[∑z-∑z*βv,∑z*βv],[vl,vv];sort = false)
    end
end

function qflash_pure(model,spec::F,x,βv,z) where F
    if spec == pressure
        T,vl,vv = saturation_temperature(model,x)
        p = x
    elseif spec == temperature
        p,vl,vv = saturation_pressure(model,x)
        T = x
    else
        throw(ArgumentError("invalid specification for qflash_pure: $spec"))
    end
    ∑z = sum(z)
    x1 = SA[1.0]

    #over critical point, or bad input
    if !isfinite(βv) || !isfinite(p) || !isfinite(T)
        return FlashResultInvalid(x1,βv)
    elseif isone(primalval(βv))
        return FlashResult([x1],[∑z*oneunit(vv)],[vv],FlashData(p,T))
    elseif iszero(primalval(βv))
        return FlashResult([x1],[∑z*oneunit(vv)],[vl],FlashData(p,T))
    elseif βv < 0 || βv > 1
        throw(error("invalid specification of vapour fraction, it must be between 0 and 1."))
    else
        return FlashResult(model,p,T,[x1,x1],[∑z-∑z*βv,∑z*βv],[vl,vv];sort = false)
    end
end

export GeneralizedXYFlash,xy_flash
