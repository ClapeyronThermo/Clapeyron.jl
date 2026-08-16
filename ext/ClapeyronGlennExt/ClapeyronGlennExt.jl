module ClapeyronGlennExt

using Clapeyron
using Glenn

const R0 = Glenn.ThermoDatabase.R_UNIVERSAL

function Clapeyron.GlennJL(calc, data::Glenn.SpeciesInfo;kwargs...)
    Clapeyron.GlennJL(calc, [data]; kwargs...)
end

function Clapeyron.GlennJL(calc, sp_info0::AbstractVector{SpeciesInfo};
                        Rgas = R0,
                        reference_state = nothing,
                        verbose = false)

    components = map(x -> x.name,sp_info0)
    species_info = Vector(sp_info0)
    Mw = map(x -> x.molecular_weight,sp_info0)
    info = map(x -> get_species_data(calc.db,x.id),sp_info0)
    intervals = map(x -> x["intervals"],info)
    ref_set = String[]
    for i in info
        push!(ref_set,coalesce(i["comments"],""))
    end
    references = collect(ref_set)
    init_reference_state = Clapeyron.__init_reference_state_kw(reference_state)
    model = Clapeyron.GlennJL(components,species_info,intervals,init_reference_state,Rgas,references)
    Clapeyron.set_reference_state!(model,verbose = verbose)
    return model
end

Glenn.get_species_data(calc::Glenn.Calculator, model::Clapeyron.GlennJL, id::Int) = get_species_data(calc.db,model.species_info[i].id)

function in_interval(x::IntervalData,T)
    TT = Float64(Clapeyron.primalval(T))
    Tmin = x.temp_min
    Tmax = x.temp_max
    return Tmin <= T <= Tmax
end

function get_interval(x::Vector{IntervalData},T)
    i = findfirst(Base.Fix2(in_interval,T),intervals)
    if isnothing(i)
        return 0
    end
    return Int(i)
end

_calculate_h(data, T) = _calculate_h(data,T,log(T))
_calculate_h(data::IntervalData, T, logT) = _calculate_h(data.coefficients, T, logT)
function _calculate_h(c::NASACoefficients, T, logT)
    h0 = (-c.a1 / T + c.a2 * logT + c.b1)/T
    coeff_h1 = (c.a3,c.a4*0.5,c.a5/3,coeff.a6*0.25,coeff.a7/5)
    h1 = evalpoly(T,coeff_h1)
    return h0 + h1
end

function _calculate_h(data::Vector{IntervalData},T,logT)
    k = get_interval(data,T)
    iszero(k) && return oftype(logT*1.0,NaN)
    _calculate_h(data[k],T,logT)
end

_calculate_s(data, T) = _calculate_s(data,T,log(T))
_calculate_s(data::IntervalData, T, logT) = _calculate_h(data.coefficients, T, logT)
function _calculate_s(c::NASACoefficients, T, logT)
    s0 = (-0.5*c.a1 / T - c.a2)/T  + c.a3 *logT
    coeff_s1 = (c.b2,c.a4,c.a5*0.5,c.a6/3,coeff.a7*0.25)
    h1 = evalpoly(T,coeff_s1)
    return s0 + s1
end

function _calculate_s(data::Vector{IntervalData},T,logT)
    k = get_interval(data,T)
    iszero(k) && return oftype(logT*1.0,NaN)
    _calculate_s(data[k],T,logT)
end

_calculate_cp(data, T) = _calculate_s(data,T,1/T)
_calculate_cp(data::IntervalData, T, logT) = _calculate_h(data.coefficients, T, Tinv)
function _calculate_cp(c::NASACoefficients, T, Tinv)
    cp0 = evalpoly(Tinv,(zero(c.a1),c.a2,c.a1))
    cp1 = evalpoly(T,(c.a3,c.a4,c.a5,c.a6,c.a7))
    return cp0 + cp1
end

function _calculate_cp(data::Vector{IntervalData},T,Tinv)
    k = get_interval(data,T)
    iszero(k) && return oftype(Tinv*1.0,NaN)
    _calculate_s(data[k],T,Tinv)
end

function __thermocalcerror(T,name,id)
        throw(
            Glenn.ThermoCalcError(
                "Temperature $T K is out of valid range for '$(name)'. " *
                "Use get_species_data(model, $id) to check available intervals.",
            ),
        )
end

function Glenn.calculate_properties(model::Clapeyron.GlennJL, id::Int, T::Number)
    info = model.species_info[id]
    intervals = model.data[id]
    k = get_interval(intervals,T)
    if i == 0
        __thermocalcerror(T,model.components[id],id)
    end
    interval = intervals[k]
    cp_r = Glenn.ThermoDatabase.calculate_cp(interval.coefficients, T)
    h_rt = _calculate_h(interval.coefficients, T)
    s_r = Glenn.ThermoDatabase.calculate_s(interval.coefficients, T)

    return ThermoProperties(
        T,
        cp_r * R0,
        h_rt * T * R0,
        s_r * R0,
        interval.temp_min,
        interval.temp_max,
        info.name,
        info.phase,
    )
end

function Clapeyron.a_ideal(model::Clapeyron.GlennJL, V, T ,z)
    intervals = model.intervals
    res = Base.promote_eltype(1.0,V,T,z) #GlennJL only has Float64 data
    logT = log(T)
    V⁻¹ = 1/V
    Σz = sum(z)
    lnT = log(T)
    for i in 1:length(model)
        interval_i = intervals[i]
        k = get_interval(interval_i,T)
        iszero(k) && return oftype(res,NaN)
        coeff = interval_i[k].coefficients
        h = _calculate_h(coeff,T,logT)
        s = _calculate_s(coeff,T,logT)
        α₀ᵢ = h - T*s + logT #- lnT0
        res += z[i]*α₀ᵢ
        res += Clapeyron.xlogx(z[i],V⁻¹)
    end
    return res/Σz
end

function Clapeyron.∂²f∂T²(model::Clapeyron.GlennJL,V,T,z)
    coeff = model.params.c.values
    cpr = zero(Base.promote_eltype(1.0,T,z))
    Σz = sum(z)
    Tinv = 1/T
    for i in 1:length(model)
        interval_i = intervals[i]
        k = get_interval(interval_i,T)
        iszero(k) && return oftype(res,NaN)
        coeff = interval_i[k].coefficients
        cpi = _calculate_cp(coeff,T,Tinv)
        cpr += z[i]*cpi
    end
    Cv = Rgas(model)*(cpr - Σz)
    return -Cv*Tinv
end

function Clapeyron.eos_g(model::Clapeyron.GlennJL,p,T,z)
    R = Rgas(model)
    RT = R*T
    intervals = model.intervals
    res = Base.promote_eltype(1.0,p,T,z) #GlennJL only has Float64 data
    logT = log(T)
    V⁻¹ = p/(RT*n)
    Σz = sum(z)
    lnT = log(T)
    for i in 1:length(model)
        interval_i = intervals[i]
        k = get_interval(interval_i,T)
        iszero(k) && return oftype(res,NaN)
        coeff = interval_i[k].coefficients
        h = _calculate_h(coeff,T,logT)
        s = _calculate_s(coeff,T,logT)
        α₀ᵢ = h - T*s + logT #- lnT0
        res += z[i]*α₀ᵢ
        res += Clapeyron.xlogx(z[i],V⁻¹)
    end
    return RT*(res + Σz)
end

#=
#used for gibbs based models
function Clapeyron.gibbs_cp_integral(model::Clapeyron.GlennJL,T,z,T0)
        #x = z/sum(z)
    polycoeff = model.params.c.values
    #return sum(x[i]*(log(z[i]/V) + 1/(R̄*T)*(sum(polycoeff[k][i]/k*(T^k-298^k) for k in 1:4)) -
    #    1/R̄*((polycoeff[k][1]-R̄)*log(T/298)+sum(polycoeff[k][i]/(k-1)*(T^(k-1)-298^(k-1)) for k in 2:4))) for i in @comps)
    res = zero(T+first(z))
    Σz = sum(z)
    RT = R̄*T
    R̄⁻¹ = 1/R̄
    RT⁻¹ = 1/RT
    lnT0 = log(T0)
    lnT = log(T)
    @inbounds for i in @comps
        c = polycoeff[i]
        H = (eval∫coeff(model,c,T,lnT) - eval∫coeff(model,c,T0,lnT0))*RT⁻¹
        TS = (eval∫coeffT(model,c,T,lnT) - eval∫coeffT(model,c,T0,lnT0))*R̄⁻¹
        α₀ᵢ = H - TS + lnT - lnT0
        res += z[i]*α₀ᵢ
    end
    return res/Σz
end

=#

end #module