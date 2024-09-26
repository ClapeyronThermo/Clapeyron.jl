struct EoSPhase{𝕋,𝕍<:Union{Nothing,AbstractVector}}
    v::𝕋 #molar volume
    t::𝕋 #temperature
    a::𝕋 #molar helmholtz energy
    p::𝕋 #pressure
    s::𝕋 #molar entropy
    sumx::𝕋 #amount of moles
    x::𝕍 #molar fractions
    phase::Symbol
    phase_type::Symbol
end

#three options for phase type
#:defined: the phase is fully defined (bulk calculations, flash)
#:zero: it is assumed that the values are molar values, but zero amount (bubble_pressure)
#:unknown: we dont know the amounts. (saturation_pressure,VLLE calculation)
#:ignore: used as placeholder for inexistent bulk phases.
function EoSPhase(model::EoSModel,V,T,z;phase = :unknown,phase_type = :defined)
    ∑z = sum(z)
    A, ∂A∂V, ∂A∂T = ∂f_vec(model,V,T,z)
    v,t,_ = promote(V/∑z,T,A)
    x = z ./ sum(z) .* one(t)
    p = -∂A∂V
    s = -∂A∂T/∑z
    a = A/∑z
    ∑x = ∑z*one(p)
    return EoSPhase(v,t,a,p,s,∑x,x,phase,phase_type)
end

EoSPhase() = EoSPhase{Nothing,Nothing}(nothing,nothing,nothing,nothing,nothing,nothing,nothing,:unknown,:ignore)


function _amount_modifier(phase,val)
    if phase.phase_type == :defined
        return val*phase.sumx
    else
        return val
    end
end

#accessors
temperature(phase::EoSPhase) = phase.t
pressure(phase::EoSPhase) = phase.p
entropy(phase::EoSPhase) = _amount_modifier(phase,phase.s)
volume(phase::EoSPhase) = _amount_modifier(phase,phase.v)
helmholtz_free_energy(phase::EoSPhase) = _amount_modifier(phase,phase.a)

function enthalpy(phase::EoSPhase)
    a,v,t,s,p = phase.a,phase.v,phase.t,phase.s,phase.p
    return _amount_modifier(phase,a + p*v + t*s)
end

function gibbs_free_energy(phase::EoSPhase)
    a,v,p = phase.a,phase.v,phase.p
    return _amount_modifier(phase,a + p*v)
end

function internal_energy(phase::EoSPhase)
    a,t,s = phase.a,phase.t,phase.s
    return _amount_modifier(phase,a + t*s)
end

function Base.show(io::IO, phase::EoSPhase)
    t,p,x = phase.t,phase.p,phase.x
    _phase = phase.phase
    if phase.phase_type != :ignore
        print(io,"EoSPhase(p = $p, t = $t, x = $x, phase = $_phase)")
    else
        print(io,"EoSPhase()")
    end
end

function Base.show(io::IO, ::MIME"text/plain", phase::EoSPhase)
    t,p,x = phase.t,phase.p,phase.x
    _phase = phase.phase
    print(io,(typeof(phase)))
    if phase.phase_type == :ignore
        print(io,"()")
        return nothing
    end
    println(io," with properties:")
    separator = " => "
    keys = (:p,:t,:x,:phase)
    vals = getfield.(Ref(phase),keys)
    #vals = [ifelse(m,missing,v) for (m,v) in zip(param.ismissingvalues, param.values)]
    show_pairs(io,keys,vals,separator,quote_string = false)
    #print(io,"EoSPhase(p = $p, t = $t, x = $x, phase = $_phase)")
end

struct EoSPhases{𝕋b,𝕍b,𝕋,𝕍,𝕍𝕍 <: AbstractVector{EoSPhase{𝕋,𝕍}}}
    bulk::EoSPhase{𝕋b,𝕍b} #bulk phase, if specified
    phases::𝕍𝕍 #vector of phases
end

Base.getindex(phases::EoSPhases,i) = phases.phases[i]
Base.length(phases::EoSPhases) = length(phases.phases)
Base.size(phases::EoSPhases) = size(phases.phases)

function EoSPhases(phases::VV) where VV <: AbstractVector{T} where T <: EoSPhase
    return EoSPhases(EoSPhase(),phases)
end

function Base.showarg(io::IO, arr::EoSPhases, toplevel)
    print(io,"EoSPhases(")
    Base.showarg(io, arr.phases, false)
    print(io, ')')
    toplevel && print(io, " with $(length(arr)) phases:")
end

function Base.show(io::IO, ::MIME"text/plain", phases::EoSPhases)
    Base.showarg(io::IO, phases, true)
    println(io)
    Base.print_matrix(IOContext(io, :compact => true),phases.phases)
end

export EoSPhase, EoSPhases
#=

z^(x/y) = z^(x)^(1/y) 
=#