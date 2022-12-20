import ThermoState
import ThermoState: @spec_str, FromState
import Unitful: @u_str

Unitful.@derived_dimension __MassDensity Unitful.𝐌/Unitful.𝐋^3
Unitful.@derived_dimension __MolDensity Unitful.𝐍/Unitful.𝐋^3
Unitful.@derived_dimension __MassVolume Unitful.𝐋^3/Unitful.𝐌
Unitful.@derived_dimension __MolVolume Unitful.𝐋^3/Unitful.𝐍

const __VolumeKind = Union{__MassDensity,__MolDensity,__MassVolume,__MolVolume,Unitful.Volume}

#real, use st argument
autospec(λ::Real,st) = ThermoState.spec(st,λ)



#with units, use units
#p,T
autospec(λ::Unitful.Temperature,st) = ThermoState.spec(spec"T",λ)
autospec(λ::Unitful.Pressure,st) = ThermoState.spec(spec"p",λ)
autospec(λ::Unitful.Volume,st) = ThermoState.spec(spec"total_v",λ)
autospec(λ::__MolDensity,st) = ThermoState.spec(spec"mol_rho",λ)
autospec(λ::__MassDensity,st) = hermoState.spec(spec"mass_rho",λ)
autospec(λ::__MolVolume,st) = ThermoState.spec(spec"mol_v",λ)
autospec(λ::__MassVolume,st) = ThermoState.spec(spec"mass_v",λ)

#mass, moles

#in case of passing only a number instead of a vector in z
autospec(t::Unitful.Mass,st) = ThermoState.spec(spec"mass",t)
autospec(t::Unitful.Amount,st) = ThermoState.spec(spec"moles",t)

#Vectors of mass, moles

autospec(t::AbstractVector{<:Unitful.Amount},st) = isone(length(t)) ? ThermoState.spec(spec"moles",only(t)) : ThermoState.spec(spec"n",t)
autospec(t::AbstractVector{<:Unitful.Mass},st) = isone(length(t)) ? ThermoState.spec(spec"mass",only(t)) : ThermoState.spec(spec"m",t)

#Clapeyron default: if no numbers are passed, assume molar amounts
autospec(t::AbstractVector{<:Real},st) = isone(length(t)) ? ThermoState.spec(spec"moles",only(t)) : ThermoState.spec(spec"n",t)

ThermoState.molecular_weight(model::EoSModel) = mw(model)
#transforms vt or pt model to ThermoState.state
function standarize(model,λ,t,z)
    spec_t = autospec(t,spec"T")
    spec_λ = autospec(λ,spec"p") #using pressure by default
    spec_z = autospec(z,spec"n")
    st = ThermoState.state(spec_λ,spec_t,spec_z)
end

#state_to_vt:  converts a ThermoState.state to V,T,z valid Clapeyron coordinates
state_to_vt(model,st::ThermoState.ThermodynamicState) = state_to_vt(ThermoState.state_type(st),model,st)

function state_to_vt(::ThermoState.QuickStates.SingleVT,model,st::ThermoState.ThermodynamicState)
    mww = only(mw(model))
    v = ThermoState.total_volume(FromState(),st,u"m^3",mww) 
    T = ThermoState.temperature(FromState(),st,u"K")
    _z = ThermoState.moles(FromState(),st,u"mol",mww)
    #_z = ThermoState.moles(FromState(),st,mw(model))
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

function state_to_pt(::ThermoState.QuickStates.SinglePT,model,st::ThermoState.ThermodynamicState)
    p = ThermoState.pressure(FromState(),st,u"Pa") |> only
    T = ThermoState.temperature(FromState(),st,u"K") |> only
    _z = ThermoState.moles(FromState(),st,u"mol",mw(model))
    z = SA[_z]
    return p,T,z
end
function state_to_pt(::ThermoState.QuickStates.MultiPT,model,st::ThermoState.ThermodynamicState)
    p = ThermoState.pressure(FromState(),st,u"Pa") |> only
    T = ThermoState.temperature(FromState(),st,u"K") |> only
    z = ThermoState.mol_number(FromState(),st,u"mol",mw(model))
    return p,T,z
end



