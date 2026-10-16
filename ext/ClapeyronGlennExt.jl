module ClapeyronGlennExt

using Clapeyron
using Glenn

const R0 = Glenn.ThermoDatabase.R_UNIVERSAL

_vec(x::AbstractVector) = x
_vec(x) = [x]
tdb(calc::Glenn.Calculator) = calc.db
tdb(db::Glenn.ThermoDB) = db

function glenn_sp_info_from_input(calc::Glenn.Calculator,input)
    return glenn_sp_info_from_input(calc.db,_vec(input))
end

function glenn_sp_info_from_input(calc::Glenn.ThermoDB,input)
    nc = length(input)
    sp_info = Vector{SpeciesInfo}(undef,nc)


    if input isa AbstractVector{SpeciesInfo}
        sp_info .= input
    elseif input isa AbstractVector && eltype(input) <: AbstractString
        for i in 1:nc
            sp_info[i] = only(Glenn.find_species(calc, String(input[i]),exact_match = true))
        end
    elseif input isa AbstractVector && eltype(input) <: Integer
        for i in 1:nc
            sp_info[i] = get_species_info(calc,Int(input[i]))
        end
    else
        throw(error("GlennJL: invalid input type :$(typeof(input))"))
    end
    return sp_info
end

function Clapeyron.__glenn_jl(calc, input;
                        Rgas = Clapeyron.Rgas(),
                        R0 = calc.R_ref,
                        reference_state = nothing,
                        verbose = false,
                        strict = true)

    species_info = glenn_sp_info_from_input(calc,input)
    nc = length(species_info)
    components = map(x -> x.name,species_info)
    info = map(x -> get_species_data(tdb(calc),x.id),species_info)

    _R0 = Vector{Float64}(undef,nc)
    if R0 isa Number || R0 isa AbstractVector
        _R0 .= R0 ./ Rgas
    else
        _R0 .= 1.0
    end

    intervals = map(x -> x["intervals"],info)

    if strict && nc > 1
        state1 = species_info[1].phase
        Tmin = intervals[1][begin].temp_min
        Tmax = intervals[1][end].temp_max
        for i in 2:nc
            if state1 != species_info[i].phase
                throw(error("GlennJL: model being instantiated with different phases. This check can be ommited by passing the keyword `strict = false`."))
            end
            Tmin = max(Tmin,intervals[i][begin].temp_min)
            Tmax = min(Tmax,intervals[i][begin].temp_max)
            if Tmin >= Tmax
                throw(error("GlennJL: null joint temperature interval. This check can be ommited by passing the keyword `strict = false`."))
            end
        end
    end


    ref_set = String[]
    for i in info
        push!(ref_set,coalesce(i["comments"],""))
    end
    references = collect(ref_set)
    init_reference_state = Clapeyron.__init_reference_state_kw(reference_state)
    model = Clapeyron.GlennJL(components,species_info,intervals,init_reference_state,Rgas,_R0,references)
    Clapeyron.set_reference_state!(model,verbose = verbose)
    return model
end

function Clapeyron.molecular_weight(model::Clapeyron.GlennJL,z)
    Clapeyron.check_arraysize(model,z)
    info = model.species_info
    res = zero(Base.promote_eltype(1.0,z))
    for i in eachindex(z)
        mw = info[i].molecular_weight
        res += mw*z[i]
    end
    return res*0.001
end

function Clapeyron.mw(model::Clapeyron.GlennJL)
    info = model.species_info
    return map(x -> x.molecular_weight,info)
end

function in_interval(x::IntervalData,T)
    TT = Float64(Clapeyron.primalval(T))
    Tmin = x.temp_min
    Tmax = x.temp_max
    return Tmin <= TT <= Tmax
end

function get_interval(intervals::Vector{IntervalData},T)
    i = findfirst(Base.Fix2(in_interval,T),intervals)
    if isnothing(i)
        return 0
    end
    return Int(i)
end

_calculate_h(data, T) = _calculate_h(data,T,log(T))
_calculate_h(data::IntervalData, T, logT) = _calculate_h(data.coefficients, T, logT)
function _calculate_h(c::NASACoefficients, T, logT)
    h0 = (-c.a1/T + (c.a2 * logT) + c.b1)/T
    coeff_h1 = (c.a3,c.a4*0.5,c.a5/3,c.a6*0.25,c.a7/5)
    h1 = evalpoly(T,coeff_h1)
    return h0 + h1
end

function _calculate_h(data::Vector{IntervalData},T,logT)
    k = get_interval(data,T)
    iszero(k) && return oftype(logT*1.0,NaN)
    _calculate_h(data[k],T,logT)
end

_calculate_s(data, T) = _calculate_s(data,T,log(T))
_calculate_s(data::IntervalData, T, logT) = _calculate_s(data.coefficients, T, logT)
function _calculate_s(c::NASACoefficients, T, logT)
    s0 = (-0.5*c.a1 / T - c.a2)/T  + c.a3 *logT
    coeff_s1 = (c.b2,c.a4,c.a5*0.5,c.a6/3,c.a7*0.25)
    s1 = evalpoly(T,coeff_s1)
    return s0 + s1
end

function _calculate_s(data::Vector{IntervalData},T,logT)
    k = get_interval(data,T)
    iszero(k) && return oftype(logT*1.0,NaN)
    _calculate_s(data[k],T,logT)
end
_calculate_hms(data,T) = _calculate_hms(data,T,logT)
function _calculate_hms(data::Vector{IntervalData}, T, logT)
    k = get_interval(data, T)
    iszero(k) && return oftype(logT*1.0, NaN)
    _calculate_hms(data[k], T, logT)
end

function _calculate_hms(c::NASACoefficients, T, logT)
    invT = 1/T
    negpow = invT*(-0.5*c.a1*invT + (c.a2*(logT + 1) + c.b1)) - c.a3*logT + (c.a3 - c.b2)
    coeff_poly = (-c.a4*0.5, -c.a5/6, -c.a6/12, -c.a7/20)
    poly = T*evalpoly(T, coeff_poly)
    return negpow + poly
end

_calculate_cp(data, T) = _calculate_cp(data,T,1/T)
_calculate_cp(data::IntervalData, T, Tinv) = _calculate_cp(data.coefficients, T, Tinv)
function _calculate_cp(c::NASACoefficients, T, Tinv)
    cp0 = evalpoly(Tinv,(zero(c.a1),c.a2,c.a1))

    cp1 = evalpoly(T,(c.a3,c.a4,c.a5,c.a6,c.a7))
    return cp0 + cp1
end

function _calculate_cp(data::Vector{IntervalData},T,Tinv)
    k = get_interval(data,T)
    iszero(k) && return oftype(Tinv*1.0,NaN)
    _calculate_cp(data[k],T,Tinv)
end

function __thermocalcerror(T,name,id)
        throw(
            Glenn.ThermoCalcError(
                "Temperature $T K is out of valid range for '$(name)'. " *
                "Use get_species_data(model, $id) to check available intervals.",
            ),
        )
end

function Clapeyron.a_ideal(model::Clapeyron.GlennJL, V, T ,z)
    intervals = model.intervals
    res = zero(Base.promote_eltype(1.0,V,T,z)) #GlennJL only has Float64 data
    logT = log(T)
    P0 = 100000.0
    V⁻¹ = Rgas(model)/(V*P0)
    Σz = sum(z)
    f = model.R0
    for i in 1:length(model)
        interval_i = intervals[i]
        k = get_interval(interval_i,T)
        iszero(k) && return oftype(res,NaN)
        interval_ik = interval_i[k]
        coeff = interval_ik.coefficients
        #h_rt = _calculate_h(coeff,T,logT)
        #s_r = _calculate_s(coeff,T,logT)
        hs = _calculate_hms(coeff,T,logT)*f[i] #H/RT - S/R
        α₀ᵢ = hs
        res += z[i]*α₀ᵢ
        res += Clapeyron.xlogx(z[i],V⁻¹)
    end
    return res/Σz + logT - 1
end

function Clapeyron.∂²f∂T²(model::Clapeyron.GlennJL,V,T,z)
    cpr = zero(Base.promote_eltype(1.0,T,z))
    Σz = sum(z)
    Tinv = 1/T
    f = model.R0
    intervals = model.intervals
    for i in 1:length(model)
        interval_i = intervals[i]
        k = get_interval(interval_i,T)
        iszero(k) && return oftype(res,NaN)
        coeff = interval_i[k].coefficients
        cpi = _calculate_cp(coeff,T,Tinv)*f[i]
        cpr += z[i]*cpi
    end
    Cv = Rgas(model)*(cpr - Σz)
    return -Cv*Tinv
end

#=

Glenn.jl integration

=#

for f in (:calculate_h,:calculate_s,:calculate_cp)
    f2 = Symbol(:_,f)
    @eval begin
        function Glenn.$f(model::Clapeyron.GlennJL,T)
            Clapeyron.check_arraysize(model,Clapeyron.SA[1.0])
            return $f2(model.intervals[1],T)*model.R0[1]
        end
    end

    @eval begin
        function Glenn.$f(model::Clapeyron.GlennJL,T,z)
            Clapeyron.check_arraysize(model,z)
            r = model.R0
            res = zero(Base.promote_eltype(1.0,T,z))
            for i in eachindex(z)
                res += $f2(model.intervals[i],T)*z[i]*r[i]
            end
            return res/sum(z)
        end
    end
end

const MIXTURE = "mixture"
const GAS = "gas"

function Glenn.calculate_properties(model::GlennJL, T, z = Clapeyron.SA[1.0])
    T = Float64(Clapeyron.primalval(T))
    R = Rgas(model)
    cp = calculate_cp(model,T,z) |> Clapeyron.primalval |> Float64
    h_rt = calculate_h(model,T,z) |> Clapeyron.primalval |> Float64
    s_r = calculate_s(model,T,z) |> Clapeyron.primalval |> Float64
    if length(model) == 1
        name = model.species_info[1].name
        Tmin = model.intervals[1][begin].temp_min
        Tmax = model.intervals[1][end].temp_max
        state = model.species_info[1].phase
    else
        name = MIXTURE
        Tmin = NaN
        Tmax = NaN
        state = GAS
    end

    ThermoProperties(T,cp,h_rt*R*T,s_r*R,Tmin,Tmax,name,state)
end

function Glenn.calculate_formation_enthalpy(model::GlennJL, z = Clapeyron.SA[1.0])
    res = zero(Base.promote_eltype(1.0,z))
    info = model.species_info
    for i in eachindex(z)
        info_i = info[i]
        res += info_i.heat_of_formation_298K * z[i]
    end
    return res/sum(z)
end

function Glenn.calculate_enthalpy_change(model::GlennJL, T1, T2, z = Clapeyron.SA[1.0])
    h1 = calculate_h(model,T1,z)*T1
    h2 = calculate_h(model,T2,z)*T2
    return (h2  - h1)*Rgas(model)
end

function Glenn.get_properties_range(model::GlennJL, T, z = Clapeyron.SA[1.0])
    return map(Ti -> Glenn.calculate_properties(model,Ti,z),T)
end

end #module