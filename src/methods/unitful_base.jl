import ThermoState
import ThermoState: @spec_str, FromState
import Unitful: @u_str

Unitful.@derived_dimension __MassDensity Unitful.𝐌/Unitful.𝐋^3
Unitful.@derived_dimension __MolDensity Unitful.𝐍/Unitful.𝐋^3
Unitful.@derived_dimension __MassVolume Unitful.𝐋^3/Unitful.𝐌
Unitful.@derived_dimension __MolVolume Unitful.𝐋^3/Unitful.𝐍

const __VolumeKind = Union{__MassDensity,__MolDensity,__MassVolume,__MolVolume,Unitful.Volume}

#real, use st argument
function autospec(λ::Real,st)
    return ThermoState.spec(st,λ)
end

#with units, use units
function autospec(λ::Unitful.Temperature,st)
    return ThermoState.spec(spec"T",λ)
end

function autospec(λ::Unitful.Pressure,st)
    return ThermoState.spec(spec"p",λ)
end

function autospec(λ::Unitful.Volume,st)
    return ThermoState.spec(spec"total_v",λ)
end

function autospec(λ::__MolDensity,st)
    return ThermoState.spec(spec"mol_rho",λ)
end

function autospec(λ::__MassDensity,st)
    return ThermoState.spec(spec"mass_rho",λ)
end

function autospec(λ::__MolVolume,st)
    return ThermoState.spec(spec"mol_v",λ)
end

function autospec(λ::__MassVolume,st)
    return ThermoState.spec(spec"mass_v",λ)
end

function autospec(t::AbstractVector{<:Real},st)
    if length(t) == 1
        return ThermoState.spec(spec"moles",only(t))
    else
        return ThermoState.spec(spec"n",t)
    end
end

function autospec(t::AbstractVector{<:Unitful.Amount},st)
    if length(t) == 1
        return ThermoState.spec(spec"moles",only(t))
    else
        return ThermoState.spec(spec"n",t)
    end
end

function autospec(t::AbstractVector{<:Unitful.Mass},st)
    if length(t) == 1
        return ThermoState.spec(spec"mass",only(t))
    else
        return ThermoState.spec(spec"m",t)
    end
end

ThermoState.molecular_weight(model::EoSModel) = mw(model)

#transforms vt or pt model to ThermoState.state
function standarize(model,λ,t,z)
    spec_t = autospec(t,spec"T")
    spec_λ = autospec(λ,spec"p") #using pressure by default
    spec_z = autospec(z,spec"n")
    st = ThermoState.state(spec_λ,spec_t,spec_z)
end

function state_to_vt(model,st::ThermoState.ThermodynamicState)
    return state_to_vt(ThermoState.state_type(st),model,st)
end

function state_to_vt(::ThermoState.QuickStates.SingleVT,model,st::ThermoState.ThermodynamicState)
    v = ThermoState.total_volume(FromState(),st,u"m^3",mw(model)) |> only
    T = ThermoState.temperature(FromState(),st,u"K") |> only
    _z = ThermoState.moles(FromState(),st,u"mol",mw(model))
    z = SA[_z]
    return v,T,z
end

function state_to_vt(::ThermoState.QuickStates.MultiVT,model,st::ThermoState.ThermodynamicState)
    v = ThermoState.total_volume(FromState(),st,u"m^3",mw(model)) |> only
    T = ThermoState.temperature(FromState(),st,u"K") |> only
    z = ThermoState.mol_number(FromState(),st,u"mol",mw(model))
    return v,T,z
end
function state_to_pt(model,st::ThermoState.ThermodynamicState)
    return state_to_pt(ThermoState.state_type(st),model,st)
end

function state_to_pt(::ThermoState.QuickStates.SingleVT,model,st::ThermoState.ThermodynamicState)
    p = ThermoState.pressure(FromState(),st,u"Pa") |> only
    T = ThermoState.temperature(FromState(),st,u"K") |> only
    _z = ThermoState.moles(FromState(),st,u"mol",mw(model))
    z = SA[_z]
    return p,T,z
end
function state_to_pt(::ThermoState.QuickStates.MultiVT,model,st::ThermoState.ThermodynamicState)
    p = ThermoState.pressure(FromState(),st,u"Pa") |> only
    T = ThermoState.temperature(FromState(),st,u"K") |> only
    z = ThermoState.mol_number(FromState(),st,u"mol",mw(model))
    return p,T,z
end



