## LLE pressure solver
function __x0_LLE_pressure(model::EoSModel,T,x,p0 = nothing,cache = nothing)
    if p0 != nothing
        p0x = p0
        vx = volume(model,p0x,T,x,phase = :l)
    elseif __γ_unwrap(model) isa ActivityModel && !(__γ_unwrap(model) isa IdealLiquidSolution)
        p0x = zero(Base.promote_eltype(model,T,x))
        vx = volume(model,p0x,T,x,phase = :l)
    else
        p0x,vx,_ = edge_pressure(model,T,x)
        #vx = volume(model,p0x,T,x,phase = :l)
    end

    tpd_result = tpd2(model,p0x,T,x,cache,lle = true, strategy = :pure, break_first = true, reduced = true)
    comps = tpd_result.compositions
    phases = tpd_result.phases
    if length(comps) == 1
        w = comps[1]
        vw = tpd_result.volumes[1]
    else
        w = similar(x,typeof(T0x))
        w .= NaN
        vw = oftype(T0x,NaN)
    end
    return p0x,vx,vw,w
end

"""
    LLE_pressure(model::EoSModel, T, x; v0 = x0_LLE_pressure(model,T,x))

Calculates the Liquid-Liquid equilibrium pressure and properties at a given temperature `T`.

Returns a tuple, containing:
- LLE Pressure `[Pa]`
- Liquid molar volume of composition `x₁ = x` at LLE Point `[m³·mol⁻¹]`
- Liquid molar volume of composition `x₂` at LLE Point  `[m³·mol⁻¹]`
- Liquid composition `x₂`
"""
function LLE_pressure(model::EoSModel,T,x;kwargs...)
    if keys(kwargs) == (:v0,)
        nt_kwargs = NamedTuple(kwargs)
        v0 = nt_kwargs.v0
        vx = exp10(v0[1])
        vw = exp10(v0[2])
        vol0 = (vx,vw)
        w0 = v0[3:end]
        _kwargs = (;vol0,w0)
        method = init_preferred_method(LLE_pressure,model,_kwargs)
    else
        method = init_preferred_method(LLE_pressure,model,kwargs)
    end
    return LLE_pressure(model, T, x, method)
end

function LLE_pressure(model::EoSModel, T, x, method::ThermodynamicMethod)
    moles_positivity(x)
    x = x/sum(x)
    T = float(T)
    verbose = get_verbosity(method)
    _model_r,idx_r = index_reduction(model,x)
    multiple_component_check(method,_model_r)
    model_r = __tpflash_cache_model(_model_r,p,NaN,x,:lle)
    x_r = x[idx_r]
    method_r = index_reduction(method,idx_r)
    λmodel,λT,λx = primalval(model_r),primalval(T),primalval(x_r)
    λresult = LLE_pressure_impl(λmodel,λT,λx,primalval(method_r))
    tup = (model_r,T,x_r)
    if any(has_dual,tup)
        λtup = (λmodel,λT,λx)
        result = lle_pressure_ad(λresult,tup,λtup)
    else
        result = λresult
    end

    (P_sat, v_x, v_w, w_r) = result
    w = index_expansion(w_r,idx_r)
    converged = bubbledew_check(model,P_sat,T,v_w,v_x,w,x)

verbose && @info "LLE_pressure results:
p  = $(primalval(P_sat))
vx = $(primalval(v_x))
vw = $(primalval(v_w))
w  = $(primalval(w))"

verbose && !converged && @info "LLE_pressure: convergence checks failed."
    if converged
        return (P_sat, v_x, v_w, w)
    else
        nan = zero(v_w)/zero(v_w)
        w = w*nan
        return (nan,nan,nan,w)
    end
end

"""
    LLE_temperature(model::EoSModel, p, x; T0 = x0_LLE_temperature(model,p,x))

Calculates the Liquid-Liquid equilibrium temperature and properties at a given pressure `p`.

Returns a tuple, containing:
- LLE Pressure `[Pa]`
- Liquid molar volume of composition `x₁ = x` at LLE Point `[m³·mol⁻¹]`
- Liquid molar volume of composition `x₂` at LLE Point  `[m³·mol⁻¹]`
- Liquid composition `x₂`
"""
function LLE_temperature(model::EoSModel,p,x;kwargs...)
    moles_positivity(x)
    if keys(kwargs) == (:v0,)
        nt_kwargs = NamedTuple(kwargs)
        v0 = nt_kwargs.v0
        T0 = v0[1]
        vl = exp10(v0[2])
        vv = exp10(v0[3])
        vol0 = (vl,vv)
        w0 = v0[4:end]
        _kwargs = (;T0,vol0,w0)
        method = init_preferred_method(LLE_temperature,model,_kwargs)
    else
        method = init_preferred_method(LLE_temperature,model,kwargs)
    end
    return LLE_temperature(model,p,x,method)
end

function LLE_temperature(model::EoSModel, p, x, method::ThermodynamicMethod)
    moles_positivity(x)
    x = x/sum(x)
    p = float(p)
    verbose = get_verbosity(method)
    _model_r,idx_r = index_reduction(model,x)
    multiple_component_check(method,_model_r)
    x_r = x[idx_r]
    model_r = __tpflash_cache_model(_model_r,p,NaN,x,:lle)
    method_r = index_reduction(method,idx_r)
    λmodel,λp,λx = primalval(model_r),primalval(p),primalval(x_r)
    λresult = LLE_temperature_impl(λmodel,λp,λx,primalval(method_r))
    tup = (model_r,p,x_r)
    if any(has_dual,tup)
        λtup = (λmodel,λp,λx)
        result = lle_temperature_ad(λresult,tup,λtup)
    else
        result = λresult
    end

    (T_sat, v_x, v_w, w_r) = result
    w = index_expansion(w_r,idx_r)
    converged = bubbledew_check(model,p,T_sat,v_w,v_x,w,x)

    verbose && @info "LLE_temperature results:
    T  = $(primalval(T_sat))
    vz = $(primalval(v_x))
    vw = $(primalval(v_w))
    w  = $(primalval(w))"

    verbose && !converged && @info "LLE_temperature: convergence checks failed."

    if converged
        return (T_sat, v_x, v_w, w)
    else
        nan = zero(v_w)/zero(v_w)
        y = y*nan
        return (nan,nan,nan,w)
    end
end

function __x0_LLE_temperature(model::EoSModel,p,x,T0 = nothing)
    if model isa PTFlashWrapper

        #check if we can actually use this here
    end
    
    if T0 != nothing
        T0x = T0
        vx = volume(model,p,T0x,x,phase = :l)
    else
        T0x,vx,_ = _edge_temperature(model,p,x)
    end

    tpd_result = tpd2(model,p,T0x,x,lle = true, strategy = :pure, break_first = true, reduced = true)
    comps = tpd_result.compositions
    phases = tpd_result.phases
    if length(comps) == 1
        w = comps[1]
        vw = tpd_result.volumes[1]
    else
        w = similar(x,typeof(T0x))
        w .= NaN
        vw = oftype(T0x,NaN)
    end
    return T0x,vlx,vlw,w
end

function LLE_pressure_init(model,T,x,vol0,p0,y0,volatiles = FillArrays.Fill(true,length(model)),verbose = false)
    if !isnothing(y0)
        if !isnothing(p0)
            if !isnothing(vol0)
                vl,vv = vol0
                verbose && @info "LLE_pressure: pressure,volumes and compositions already provided."
            else
                verbose && @info "LLE_pressure: calculating volumes from provided pressure and compositions."
                vl = volume(model,p0,T,x,phase = :l)
                vv = volume(model,p0,T,y0,phase = :l)
            end
        else
            if !isnothing(vol0)
                vl,vv = vol0
                verbose && @info "LLE_pressure: calculating pressure from provided vapour volume and composition."
                if model isa PTFlashWrapper
                    p0 = vl*NaN
                else
                    p0 = pressure(model,vv,T,y0)
                end
            else
                verbose && @info "LLE_pressure: calculating volumes and pressures from provided vapour composition."
                p0,_,_,_ = __x0_LLE_pressure(model,T,x)
                vl = volume(model,p0,T,x,phase = :l)
                vv = volume(model,p0,T,y0,phase = :l)
            end
        end
    else
        p00,vl0,vv0,y0 = __x0_LLE_pressure(model,T,x,p0)
        if !isnothing(p0)
            verbose && @info "LLE_pressure: calculating volumes and compositions from provided pressure"
            vl = volume(model,p0,T,x,phase = :l)
            vv = volume(model,p0,T,y0,phase = :l)
        else
            verbose && @info "LLE_pressure: temperatures, volumes and compositions calculated from Clapeyron.__x0_LLE_presure"
            vl = vl0
            vv = vv0
            p0 = p00
        end
    end
    verbose && @info "LLE_pressure initial points:
p0: $p0
vz: $vl
vw: $vv
w0: $y0"
    return p0,vl,vv,y0
end

