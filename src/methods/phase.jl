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

EosPhase() = EoSPhase(nothing,nothing,nothing,nothing,nothing,nothing,nothing,:unknown,:unknown)

temperature(phase::EoSPhase) = phase.T

function _amount_modifier(phase,val)
    if phase.phase_type == :defined
        return val*phase.sumx
    else
        return val
    end
end

#accessors
pressure(phase::EoSPhase) = phase.p
entropy(phase::EoSPhase) = _amount_modifier(phase,phase.s)
volume(phase::EoSPhase) = _amount_modifier(phase,phase.v)
helmholtz_free_energy(phase::EoSPhase) = _amount_modifier(phase,phase.a)

function enthalpy(phase::EoSPhase)
    a,v,t,s,p = phase.a,phase.v,phase.t,phase.s,phase.p
    return _amount_modifier(phase,a + p*v + t*s)
end

function gibbs_free_energy(phase::EoSPhase)
    a,v,t,s,p = phase.a,phase.v,phase.p
    return _amount_modifier(phase,a + p*v)
end

function internal_energy(phase::EoSPhase)
    a,v,t,s,p = phase.a,phase.v,phase.t,phase.s,phase.p
    return _amount_modifier(phase,a + t*s)
end

function Base.show(io::IO, phase::EoSPhase)
    t,p,x = phase.t,phase.p,phase.x
    _phase = phase.phase
    print(io,"EoSPhase(p = $p, t = $t, x = $x, phase = $_phase)")
end

function Base.show(io::IO, ::MIME"text/plain", phase::EoSPhase)
    t,p,x = phase.t,phase.p,phase.x
    _phase = phase.phase
    println(io(typeof(phase)),"with properties:")
    
    print(io,"EoSPhase(p = $p, t = $t, x = $x, phase = $_phase)")
end

function EoSPhases{𝕋b,𝕍b,𝕋,𝕍}
    bulk::EoSPhase{𝕋b,𝕍b} #bulk phase, if specified
    phases::Vector{EoSPhase{𝕋,𝕍}} #vector of phases
end

function EoSPhases(phases::Vector{EoSPhase{𝕋,𝕍}})
    return EoSPhases(EoSPhase(),phases)
end

