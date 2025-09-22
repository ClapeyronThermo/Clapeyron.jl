module ClapeyronJutulDarcyExt
    using Clapeyron: Clapeyron
    using JutulDarcy: JutulDarcy, MultiComponentFlash, Jutul
    const C = Clapeyron
    const M = MultiComponentFlash
    const J = JutulDarcy
    function J.MultiPhaseCompositionalSystemLV(
            equation_of_state::C.EoSModel, 
            phases = (J.LiquidPhase(), J.VaporPhase()); 
            reference_densities = ones(length(phases)), 
            other_name = "Water",
            reference_phase_index = J.get_reference_phase_index(phases))

        
        c = copy(equation_of_state.components)
        N = length(c)
        phases = tuple(phases...)
        T = typeof(phases)
        nph = length(phases)
        @assert nph == 2 || nph == 3
        reference_densities = tuple(reference_densities...)
        @assert length(reference_densities) == nph
        if nph == 3
            other = only(filter(x -> !(isa(x, LiquidPhase) || isa(x, VaporPhase)), phases))
            O = typeof(other)
            push!(c, other_name)
        else
            O = Nothing
        end
        only(findall(isequal(J.LiquidPhase()), phases))
        only(findall(isequal(J.VaporPhase()), phases))
        J.MultiPhaseCompositionalSystemLV{typeof(equation_of_state), T, O, typeof(reference_densities),N}(phases, c, equation_of_state, reference_densities, reference_phase_index)
    end

    function Base.show(io::IO, sys::J.MultiPhaseCompositionalSystemLV{<:C.EoSModel})
        n = length(sys.equation_of_state)
        if J.has_other_phase(sys)
            name = "(with water)"
            n = n - 1
        else
            name = "(no water)"
        end
        eos = sys.equation_of_state
        cnames = join(C.component_list(eos), ", ")
        print(io, "MultiPhaseCompositionalSystemLV $name with $(Base.summary(eos)) EOS with $n components: $cnames")
    end

    J.get_compressibility_factor(forces, eos::C.EoSModel, P, T, Z,phase = :unknown) = C.compressibility_factor(eos, P, T, Z,phase = phase)
    
    @inline function J.single_phase_update!(P, T, Z, x, y, forces, eos::C.EoSModel, c)
        AD = Base.promote_type(eltype(Z), typeof(P), typeof(T))
        @. x = Z
        @. y = Z
        V = M.single_phase_label(eos, c)
        if V > 0.5
            phase_state = MultiComponentFlash.single_phase_v
            Z_V = J.get_compressibility_factor(forces, eos, P, T, Z, :v)
            Z_L = Z_V
        else
            phase_state = MultiComponentFlash.single_phase_l
            Z_L = J.get_compressibility_factor(forces, eos, P, T, Z, :l)
            Z_V = Z_L
        end
        V = convert(AD, V)
        out = (Z_L::AD, Z_V::AD, V::AD, phase_state::M.PhaseState2Phase)
        return out
    end

    @inline function J.two_phase_update!(S, P, T, Z, x, y, K, vapor_frac, forces, eos::C.EoSModel, c)
        AD = Base.promote_type(typeof(P), eltype(Z), typeof(T))
        @. x = J.liquid_mole_fraction(Z, K, vapor_frac)
        @. y = J.vapor_mole_fraction(x, K)
        V = J.two_phase_pre!(S, P, T, Z, x, y, vapor_frac, eos, c)
        Z_L = J.get_compressibility_factor(forces, eos, P, T, x, :l)
        Z_V = J.get_compressibility_factor(forces, eos, P, T, y, :v)
        phase_state = MultiComponentFlash.two_phase_lv
    
        return (Z_L::AD, Z_V::AD, V::AD, phase_state)
    end

#=
const FD = C.ForwardDiff

FD.:≺(::J.Cells, ::Type{<:FD.Tag}) = false
FD.:≺(::Type{<:FD.Tag}, ::J.Cells) = true
=#
end #module
