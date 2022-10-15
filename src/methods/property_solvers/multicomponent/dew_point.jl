
"""
    DewPointMethod <: ThermodynamicMethod 

Abstract type for `dew_pressure` and `dew_temperature` routines.

Should at least support passing the `x0` keyword, containing an initial vapour phase, if available.

"""
abstract type DewPointMethod <: ThermodynamicMethod end


function index_reduction(method::DewPointMethod,idx_r)
    if hasfield(typeof(method),:x0)
        if !isnothing(method.x0) 
            method_r = deepcopy(method)
            x0_new = method.x0[idx_r]
            resize!(method_r.x0,length(x0_new))
            method_r.x0 .= x0_new
            return method_r
        end
    end
    return method
end

function __x0_dew_pressure(model::EoSModel,T,y)
    comps = length(model)
    pure = split_model(model)
    crit = crit_pure.(pure)
    
    T_c = [tup[1] for tup in crit]
    V_c = [tup[3] for tup in crit]
    _0 = zero(T+first(y))
    nan = _0/_0 
    sat_nan = (nan,nan,nan)
    sat = fill(sat_nan,comps)
    replaceP = T_c .< T
    for i in 1:comps
        if !replaceP[i]
        sat[i] = saturation_pressure(pure[i],T,ChemPotVSaturation(crit = crit[i]))
        end
    end
    
    P_sat = [tup[1] for tup in sat]
    V_l_sat = [tup[2] for tup in sat]
    V_v_sat = [tup[3] for tup in sat]
    P⁻¹ = zero(T)
    V0_l = zero(T)
    V0_v = zero(T)
    Pi   = zero(y)
    for i in 1:length(y)
        if !replaceP[i]
            Pi[i] = P_sat[i][1]
            P⁻¹+=y[i]/Pi[i]
            V0_v += y[i]*V_v_sat[i]
        else 
            Pi[i] = pressure(pure[i],V_c[i],T)
            P⁻¹+=y[i]/Pi[i]
            V0_v += y[i]*V_c[i]*1.2
        end
    end
    P = 1/P⁻¹
    x = @. y*P/Pi
    xsum = 1/∑(x)
    x    = x.*xsum
    
    for i in 1:length(y)
        if !replaceP[i]
            V0_l += x[i]*V_l_sat[i]
        else
            V0_l += x[i]*V_c[i]
        end
    end
    return P,V0_l,V0_v,x
    #prepend!(x,log10.([V0_l,V0_v]))
    #return x
end

function x0_dew_pressure(model::EoSModel,T,y)
    P,V0_l,V0_v,x = __x0_dew_pressure(model,T,y)
    prepend!(x,log10.([V0_l,V0_v]))
    return x
end

function dew_pressure_init(model,T,y,vol0,p0,x0)
    if !isnothing(x0)
        if !isnothing(p0)
            if !isnothing(vol0)
                vl,vv = vol0
            else
                vl = volume(model,p0,T,x0,phase = :l)
                vv = volume(model,p0,T,y,phase =:v)
            end
        else
            if !isnothing(vol0)
                vl,vv = vol0
                p0 = pressure(model,vv,T,y)
            else
                p0,vl0,vv0,_ = __x0_dew_pressure(model,T,y)
                vl = min(vl0,volume(model,p0,T,x0,phase = :l))
                vv = max(vv0,volume(model,p0,T,y,phase =:v))
            end
        end
    else
        p00,vl0,vv0,x0 = __x0_dew_pressure(model,T,y)
        if !isnothing(p0)
            vl = min(vl0,volume(model,p0,T,x0,phase = :l))
            vv = max(vv0,volume(model,p0,T,y,phase = :v))
        else
            vl = vl0
            vv = vv0
            p0 = p00
        end
    end
    return p0,vl,vv,x0
end

"""
    dew_pressure(model::EoSModel, T, y,method = ChemPotDewPressure())

Calculates the dew pressure and properties at a given temperature.
Returns a tuple, containing:
- Dew Pressure `[Pa]`
- liquid volume at Dew Point [`m³`]
- vapour volume at Dew Point [`m³`]
- Liquid composition at Dew Point

By default, uses equality of chemical potentials, via [`ChemPotDewPressure`](@ref)
"""
function dew_pressure(model::EoSModel, T, y,method::DewPointMethod)
    y = y/sum(y)
    T = float(T)
    model_r,idx_r = index_reduction(model,y)
    if length(model_r)==1
        (P_sat,v_l,v_v) = saturation_pressure(model_r,T)
        return (P_sat,v_l,v_v,y)
    end
    y_r = y[idx_r]
    (P_sat, v_l, v_v, x_r) = dew_pressure_impl(model_r,T,y_r,index_reduction(method,idx_r))
    x = index_expansion(x_r,idx_r)
    converged = bubbledew_check(v_l,v_v,y,x)
    if converged
        return (P_sat, v_l, v_v, x)
    else
        nan = zero(v_l)/zero(v_l)
        x = x*nan
        return (nan,nan,nan,x)
    end
end

function antoine_dew(pure,T,y,crit)
    pᵢ = aprox_psat.(pure,T,crit)
    p = sum(y./pᵢ)^(-1)
    x = y.*p./pᵢ
    xsum = 1/∑(x)
    x    = x.*xsum
    return p,x
end

function __x0_dew_temperature(model::EoSModel,p,y)
    comps = length(model)
    pure = split_model(model)
    crit = crit_pure.(pure)
    
    p_c = [tup[2] for tup in crit]
    V_c = [tup[3] for tup in crit]
    _0 = zero(p+first(y))
    replaceP = p_c .< p
    T_sat = fill(_0,comps)
    V_l_sat = fill(_0,comps)
    V_v_sat = fill(_0,comps)
    for i in 1:comps
        crit_i = crit[i]
        Tci,Pci,Vci = crit_i
        if !replaceP[i]
            Ti,Vli,Vvi = saturation_temperature(pure[i],p,AntoineSaturation(crit = crit_i))
        else
            
            Ti,Vli,Vvi = Tci,Vci,1.2*Vci  
        end
        T_sat[i] = Ti
        V_l_sat[i] = Vli
        V_v_sat[i] = Vvi
    end    

    V0_l = zero(p)
    V0_v = zero(p)
    if !any(replaceP) #p < min(pci), proceed with entalphy aproximation:
        dPdTsat = VT_entropy.(pure,V_v_sat,T_sat) .- VT_entropy.(pure,V_l_sat,T_sat) ./ (V_v_sat .- V_l_sat)
        x = copy(dPdTsat)
        ##initialization for T, dew form
        #= we solve the aproximate problem of finding T such as:
        p = sum(yi*pi(T))
        where pi ≈ p + dpdt(T-T0)
        for a dew specification:
        sum(yi/(1 + dpdt(T-T0)/p)) - 1 = 0
        =#
        function f0p(T)
            res = zero(T)
            for i in eachindex(dPdTsat)
                pr_i = 1 + dPdTsat[i]*(T-T_sat[i])/p
                res += y[i]/pr_i
            end
            return res - 1
        end
        fTd = Roots.ZeroProblem(f0p,dot(T_sat,y))
        T0 = Roots.solve(fTd,Roots.Order0())
        sat = saturation_pressure.(pure,T0)
        for i in 1:comps
            V_l_sat[i] = sat[i][2]
            V_v_sat[i] = sat[i][3]
        end
        x .= y ./ first.(sat)
        x ./= sum(x) 
    else
        Tb = extrema(T_sat).*(0.9,1.1)
        
        f(T) = antoine_dew(pure,T,y,crit)[1]-p
        fT = Roots.ZeroProblem(f,Tb)
        T0 = Roots.solve(fT,Roots.Order0())
        p,x = antoine_dew(pure,T0,y,crit)
    end
    V0_l = zero(p)
    V0_v = zero(p)

    for i in 1:comps
        if !replaceP[i]
            V0_v += y[i]*V_v_sat[i]
            V0_l += x[i]*V_l_sat[i]
        else
            V0_v += y[i]*V_c[i]*1.2
            V0_l += x[i]*V_c[i]
        end
    end
    return T0,V0_l,V0_v,x
end

function x0_dew_temperature(model::EoSModel,p,y)
    T0,V0_l,V0_v,x = __x0_dew_temperature(model,p,y)
    prepend!(x,log10.([V0_l,V0_v]))
    prepend!(x,T0)
    return x
end

function dew_temperature_init(model,p,y,vol0,T0,x0)
    if !isnothing(x0)
        if !isnothing(T0)
            if !isnothing(vol0)
                vl,vv = vol0
            else
                vl = volume(model,p,T0,x0,phase = :l)
                vv = volume(model,p,T0,y,phase =:v)
            end
        else
            T0,vl0,vv0,_ = __x0_dew_temperature(model,p,y)
            if !isnothing(vol0)
                vl,vv = vol0
            else
                vl = min(vl0,volume(model,p,T0,x0,phase = :l))
                vv = max(vv0,volume(model,p,T0,y,phase =:v))
            end
        end
    else
        T00,vl0,vv0,x0 = __x0_dew_temperature(model,p,y)
        if !isnothing(T0)
            vl = min(vl0,volume(model,p,T0,x0,phase = :l))
            vv = max(vv0,volume(model,p,T0,y,phase = :v))
        else
            vl = vl0
            vv = vv0
            T0 = T00
        end
    end
    return T0,vl,vv,x0
end

"""
    dew_temperature(model::EoSModel, p, y, method = ChemPotDewTemperature())

calculates the dew temperature and properties at a given pressure.
Returns a tuple, containing:
- Dew Temperature `[K]`
- liquid volume at Dew Point [`m³`]
- vapour volume at Dew Point [`m³`]
- Liquid composition at Dew Point

By default, uses equality of chemical potentials, via [`ChemPotDewTemperature`](@ref)
"""
function dew_temperature(model::EoSModel,p,y,method::DewPointMethod)
    y = y/sum(y)
    p = float(p)
    model_r,idx_r = index_reduction(model,y)
    if length(model_r)==1
        (T_sat,v_l,v_v) = saturation_temperature(model_r,p)
        return (T_sat,v_l,v_v,y)
    end
    y_r = y[idx_r]
    (T_sat, v_l, v_v, x_r) = dew_temperature_impl(model_r,p,y_r,index_reduction(method,idx_r))
    x = index_expansion(x_r,idx_r)
    converged = bubbledew_check(v_l,v_v,y,x)
    if converged
        return (T_sat, v_l, v_v, x)
    else
        nan = zero(v_l)/zero(v_l)
        x = x*nan
        return (nan,nan,nan,x)
    end
end

include("dew_point/dew_chempot.jl")
include("dew_point/dew_fugacity.jl") 
include("dew_point/dew_activity.jl") 

function dew_pressure(model::EoSModel,T,y;v0 = nothing)
    if isnothing(v0)
        return dew_pressure(model,T,y,ChemPotDewPressure())
    else
        vl = exp10(v0[1])
        vv = exp10(v0[2])
        vol0 = (vl,vv)
        x0 = v0[3:end]
        dew_pressure(model,T,y,ChemPotDewPressure(;vol0,x0))
    end
end

function dew_temperature(model::EoSModel,p,y;v0 = nothing)
    if isnothing(v0)
        return dew_temperature(model,p,y,ChemPotDewTemperature())
    else
        T0 = v0[1]
        vl = exp10(v0[2])
        vv = exp10(v0[3])
        vol0 = (vl,vv)
        x0 = v0[4:end]
        dew_temperature(model,p,y,ChemPotDewTemperature(;T0,vol0,x0))
    end
end
