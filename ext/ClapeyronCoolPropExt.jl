module ClapeyronCoolPropExt

using Clapeyron
using CoolProp
using CoolProp.Unitful
using CoolProp: libcoolprop
function Clapeyron.coolprop_handler()
    Base.Libc.Libdl.dlopen(libcoolprop;throw_error = false)
end

using Clapeyron.StaticArrays
using Clapeyron: molecular_weight

#TODO: CoolProp-Clapeyron integration?

"""
    coolprop_crit_data(components)
returns a named tuple with critical and molecular weight data extracted from CoolProp.
"""
function Clapeyron.coolprop_crit_data(components)
    comps = Clapeyron.format_components(components)
    Tc = zeros(length(comps))
    Vc = zeros(length(comps))
    acentricfactor = zeros(length(comps))
    Pc = zeros(length(comps))
    Mw = zeros(length(comps))

    for (i,comp) in pairs(comps)
        Tc[i] = PropsSI("Tcrit",comp)
        Pc[i] = PropsSI("pcrit",comp)
        Vc[i] = 1/PropsSI("rhocrit",comp)
        Mw[i] = 1000*PropsSI("molarmass",comp)
        acentricfactor[i] =PropsSI("acentric",comp)
    end
    return (;Mw,Tc,Pc,Vc,acentricfactor)

end

# Define Parameters enum (equivalent to CoolProp's parameters)
@enum Parameters begin
    INVALID_PARAMETER = 0
    igas_constant
    imolar_mass
    iacentric_factor
    irhomolar_reducing
    irhomolar_critical
    iT_reducing
    iT_critical
    irhomass_reducing
    irhomass_critical
    iP_critical
    iP_reducing
    iT_triple
    iP_triple
    iT_min
    iT_max
    iP_max
    iP_min
    idipole_moment
    iT
    iP
    iQ
    iTau
    iDelta
    iDmolar
    iHmolar
    iSmolar
    iCpmolar
    iCp0molar
    iCvmolar
    iUmolar
    iGmolar
    iHelmholtzmolar
    iHmolar_residual
    iSmolar_residual
    iGmolar_residual
    iHmolar_idealgas
    iSmolar_idealgas
    iUmolar_idealgas
    iDmass
    iHmass
    iSmass
    iCpmass
    iCp0mass
    iCvmass
    iUmass
    iGmass
    iHelmholtzmass
    iHmass_idealgas
    iSmass_idealgas
    iUmass_idealgas
    iviscosity
    iconductivity
    isurface_tension
    iPrandtl
    ispeed_sound
    iisothermal_compressibility
    iisobaric_expansion_coefficient
    iisentropic_expansion_coefficient
    ifundamental_derivative_of_gas_dynamics
    ialphar
    idalphar_dtau_constdelta
    idalphar_ddelta_consttau
    ialpha0
    idalpha0_dtau_constdelta
    idalpha0_ddelta_consttau
    id2alpha0_ddelta2_consttau
    id3alpha0_ddelta3_consttau
    iBvirial
    iCvirial
    idBvirial_dT
    idCvirial_dT
    iZ
    iPIP
    ifraction_min
    ifraction_max
    iT_freeze
    iGWP20
    iGWP100
    iGWP500
    iFH
    iHH
    iPH
    iODP
    iPhase
    iundefined_parameter
end

function output_is_ideal(param::Parameters)
    return param in (iHmolar_idealgas,iSmolar_idealgas,iUmolar_idealgas,iHmass_idealgas,iSmass_idealgas,iUmass_idealgas)
end

# Define InputPairs enum (equivalent to CoolProp's input_pairs)

@enum InputPairs begin
    INPUT_PAIR_INVALID = 0
    QT_INPUTS
    PQ_INPUTS
    QSmolar_INPUTS
    QSmass_INPUTS
    HmolarQ_INPUTS
    HmassQ_INPUTS
    DmolarQ_INPUTS
    DmassQ_INPUTS
    PT_INPUTS
    DmassT_INPUTS
    DmolarT_INPUTS
    HmolarT_INPUTS
    HmassT_INPUTS
    SmolarT_INPUTS
    SmassT_INPUTS
    TUmolar_INPUTS
    TUmass_INPUTS
    DmassP_INPUTS
    DmolarP_INPUTS
    HmassP_INPUTS
    HmolarP_INPUTS
    PSmass_INPUTS
    PSmolar_INPUTS
    PUmass_INPUTS
    PUmolar_INPUTS
    HmassSmass_INPUTS
    HmolarSmolar_INPUTS
    SmassUmass_INPUTS
    SmolarUmolar_INPUTS
    DmassHmass_INPUTS
    DmolarHmolar_INPUTS
    DmassSmass_INPUTS
    DmolarSmolar_INPUTS
    DmassUmass_INPUTS
    DmolarUmolar_INPUTS
end


# Helper function to check if parameters match a specific pair (in any order)
function match_pair(key1::Parameters, key2::Parameters, x1::Parameters, x2::Parameters)
    matched = (key1 == x1 && key2 == x2) || (key1 == x2 && key2 == x1)
    swap = matched && (key1 != x1)  # Determine if values need swapping
    return matched, swap
end

function generate_update_pair(key1::Parameters,value1::T1,key2::Parameters,value2::T2) where {T1,T2}
    v1,v2 = promote(value1,value2)
    generate_update_pair(key1,v1,key2,v2)
end

# Main function to generate input pair and ensure consistent parameter order
function generate_update_pair(key1::Parameters, value1::T, key2::Parameters, value2::T) where T
    # Check all possible input pair combinations
    matched, swap = match_pair(key1, key2, iQ, iT)
    matched && return (QT_INPUTS, swap ? (value2, value1) : (value1, value2))

    matched, swap = match_pair(key1, key2, iP, iQ)
    matched && return (PQ_INPUTS, swap ? (value2, value1) : (value1, value2))

    matched, swap = match_pair(key1, key2, iP, iT)
    matched && return (PT_INPUTS, swap ? (value2, value1) : (value1, value2))

    matched, swap = match_pair(key1, key2, iDmolar, iT)
    matched && return (DmolarT_INPUTS, swap ? (value2, value1) : (value1, value2))

    matched, swap = match_pair(key1, key2, iHmolar, iT)
    matched && return (HmolarT_INPUTS, swap ? (value2, value1) : (value1, value2))

    matched, swap = match_pair(key1, key2, iSmolar, iT)
    matched && return (SmolarT_INPUTS, swap ? (value2, value1) : (value1, value2))

    matched, swap = match_pair(key1, key2, iT, iUmolar)
    matched && return (TUmolar_INPUTS, swap ? (value2, value1) : (value1, value2))

    matched, swap = match_pair(key1, key2, iDmolar, iHmolar)
    matched && return (DmolarHmolar_INPUTS, swap ? (value2, value1) : (value1, value2))

    matched, swap = match_pair(key1, key2, iDmolar, iSmolar)
    matched && return (DmolarSmolar_INPUTS, swap ? (value2, value1) : (value1, value2))

    matched, swap = match_pair(key1, key2, iDmolar, iUmolar)
    matched && return (DmolarUmolar_INPUTS, swap ? (value2, value1) : (value1, value2))

    matched, swap = match_pair(key1, key2, iDmolar, iP)
    matched && return (DmolarP_INPUTS, swap ? (value2, value1) : (value1, value2))

    matched, swap = match_pair(key1, key2, iDmolar, iQ)
    matched && return (DmolarQ_INPUTS, swap ? (value2, value1) : (value1, value2))

    matched, swap = match_pair(key1, key2, iHmolar, iP)
    matched && return (HmolarP_INPUTS, swap ? (value2, value1) : (value1, value2))

    matched, swap = match_pair(key1, key2, iP, iSmolar)
    matched && return (PSmolar_INPUTS, swap ? (value2, value1) : (value1, value2))

    matched, swap = match_pair(key1, key2, iP, iUmolar)
    matched && return (PUmolar_INPUTS, swap ? (value2, value1) : (value1, value2))

    matched, swap = match_pair(key1, key2, iHmolar, iSmolar)
    matched && return (HmolarSmolar_INPUTS, swap ? (value2, value1) : (value1, value2))

    matched, swap = match_pair(key1, key2, iSmolar, iUmolar)
    matched && return (SmolarUmolar_INPUTS, swap ? (value2, value1) : (value1, value2))

    # Return invalid pair if no match found
    return (INPUT_PAIR_INVALID, (value1, value2))
end

function standarize_value_and_itype(model,z,value,itype)
    mw = molecular_weight(model,z)/sum(z)
    itype == iDmass && return value/mw,iDmolar
    itype == iHmass && return value*mw,iHmolar
    itype == iUmass && return value*mw,iUmolar
    itype == iSmass && return value*mw,iSmolar
    itype == iHelmholtzmass && return value*mw,iHelmholtzmolar
    return value,itype
end

function CoolProp.PropsSI(output::AbstractString, name1::AbstractString, value1::Real, name2::AbstractString, value2::Real, fluid::EoSModel)
    return ClapeyronPropsSI(output,name1,value1,name2,value2,(fluid,Clapeyron.SA[1.0]))
end

function CoolProp.PropsSI(output::AbstractString, fluid::EoSModel)
    return ClapeyronPropsSI(output,(fluid,Clapeyron.SA[1.0]))
end

function CoolProp.PropsSI(output::AbstractString, fluid::Tuple{EoSModel,Any})
    return ClapeyronPropsSI(output,fluid)
end

function CoolProp.PropsSI(output::AbstractString, name1::AbstractString, value1::Union{Unitful.Quantity,Real}, name2::AbstractString, value2::Union{Unitful.Quantity,Real}, fluid::EoSModel)
    return ClapeyronPropsSI(output,name1,value1,name2,value2,(fluid,Clapeyron.SA[1.0]))
end

function CoolProp.PropsSI(output::AbstractString, name1::AbstractString, value1::Real, name2::AbstractString, value2::Real, fluid::Tuple{EoSModel,Any})
    return ClapeyronPropsSI(output,name1,value1,name2,value2,fluid)
end

function CoolProp.PropsSI(output::AbstractString, name1::AbstractString, value1::Union{Unitful.Quantity,Real}, name2::AbstractString, value2::Union{Unitful.Quantity,Real}, fluid::Tuple{EoSModel,Any})
    unit1 = CoolProp._get_unit(name1,false)
    unit2 = CoolProp._get_unit(name2,false)
    outunit = CoolProp._get_unit(output,false)
    return ClapeyronPropsSI(output, name1, CoolProp._si_value(unit1,value1), name2, CoolProp._si_value(unit2,value2), fluid)*outunit
end

function ClapeyronPropsSI(output::AbstractString, name1::AbstractString, value1::Real, name2::AbstractString, value2::Real, fluid::Tuple{EoSModel,Any})
    model,z = fluid
    i1 = Parameters(CoolProp.get_param_index(name1))
    i2 = Parameters(CoolProp.get_param_index(name2))
    if i1 == i2
        throw(error("cannot set two inputs of the same type"))
    end
    io = Parameters(CoolProp.get_param_index(output))

    if io == i1
        return value1
    end

    if io == i2
        return value2
    end

    value1_std,j1 = standarize_value_and_itype(model,z,value1,i1)
    value2_std,j2 = standarize_value_and_itype(model,z,value2,i2)
    res = generate_update_pair(j1,value1_std,j2,value2_std)
    itype,(x,y) = res
    
    flash = flash_by_input(model,x,y,z,itype)
    if output_is_ideal(io) && Clapeyron.numphases(flash) != 1
        throw(Clapeyron.invalid_property_multiphase_error(io,Clapeyron.numphases(flash),Clapeyron.pressure(flash),Clapeyron.temperature(flash)))
    end
    flash.fractions ./= sum(flash.fractions) #we move to a 1-mol basis
    fracs = z ./ sum(z)
    return eval_property(model,flash,fracs,io)
end

function ClapeyronPropsSI(output::AbstractString, fluid::Tuple{EoSModel,Any})
    model,z = fluid
    io = Parameters(CoolProp.get_param_index(output))
    fracs = z ./ sum(z)
    return eval_property(model,nothing,fracs,io)
end

function flash_not_implemented_error(input)
    str = string(Symbol(input))
    throw(error("Flash prodedure not implemented for $str"))
end

function property_not_implemented_error(input)
    str = string(Symbol(input))
    throw(error("Property evaluation of $str is not implemented."))
end

_vol(model,v,z) = 1/v

_molar(model,x,z) = x*Clapeyron.molecular_weight(model,z)
function _dmass(model,rhomass,z)
    molar_weight = molecular_weight(model,z)
    molar_weight/rhomass
end

function _dmolar(model,rhomolar,z)
    sum(z)/rhomolar
end

function _emolar(model,e,z)
    sum(z)*e
end

function _emass(model,e,z)
    molar_weight = molecular_weight(model,z)
    return e*molar_weight
end

function flash_by_input(model,x,y,z,itype::InputPairs)
    QT_INPUTS == itype && return Clapeyron.QT.flash(model,x,y,z)
    PQ_INPUTS == itype && return Clapeyron.QP.flash(model,y,x,z)
    QSmolar_INPUTS == itype && flash_not_implemented_error(itype)
    HmolarQ_INPUTS == itype && flash_not_implemented_error(itype)
    DmolarQ_INPUTS == itype && flash_not_implemented_error(itype)
    PT_INPUTS == itype && return Clapeyron.PT.flash(model,x,y,z)
    DmolarT_INPUTS == itype && return Clapeyron.PT.flash(_dmolar(model,x,z),y,z)
    HmolarT_INPUTS == itype && flash_not_implemented_error(itype)
    SmolarT_INPUTS == itype && return Clapeyron.TS.flash(model,y,_emolar(model,x,z),z)
    TUmolar_INPUTS == itype && flash_not_implemented_error(itype)
    DmolarP_INPUTS == itype && flash_not_implemented_error(itype)
    HmolarP_INPUTS == itype && return Clapeyron.PH.flash(model,y,_emolar(model,x,z),z)
    PSmolar_INPUTS == itype && return Clapeyron.PS.flash(model,x,_emolar(model,y,z),z)
    PUmolar_INPUTS == itype && flash_not_implemented_error(itype)
    HmolarSmolar_INPUTS == itype && flash_not_implemented_error(itype)
    SmolarUmolar_INPUTS == itype && flash_not_implemented_error(itype)
    DmolarHmolar_INPUTS == itype && flash_not_implemented_error(itype)
    DmolarSmolar_INPUTS == itype && flash_not_implemented_error(itype)
    DmolarUmolar_INPUTS == itype && Clapeyron.uv_flash(model,_emolar(model,y,z),_dmolar(model,x,z),z)
end

function eval_property(model,result,z,key::Parameters)
    igas_constant == key && return Clapeyron.Rgas(model)
    imolar_mass == key && return Clapeyron.molecular_weight(model,z)
    iacentric_factor == key && return Clapeyron.acentric_factor(model)
    irhomolar_reducing == key && return property_not_implemented_error(key)
    irhomolar_critical == key && return sum(z)/Clapeyron.crit_pure(model)[3]
    iT_reducing == key && return property_not_implemented_error(key)
    iT_critical == key && crit_pure(model)[1]
    irhomass_reducing == key && return property_not_implemented_error(key)
    irhomass_critical == key && return molecular_weight(model,SA[1.0])/crit_pure(model)[3]
    iP_critical == key && return crit_pure(model)[2]
    iP_reducing == key && return property_not_implemented_error(key)
    iT_triple == key && return property_not_implemented_error(key)
    iP_triple == key && return property_not_implemented_error(key)
    iT_min == key && return property_not_implemented_error(key)
    iT_max == key && return property_not_implemented_error(key)
    iP_max == key && return property_not_implemented_error(key)
    iP_min == key && return property_not_implemented_error(key)
    idipole_moment == key && return property_not_implemented_error(key)
    iT == key && return Clapeyron.temperature(result)
    iP == key && return Clapeyron.pressure(result)
    iQ == key && return property_not_implemented_error(key)
    iTau == key && return property_not_implemented_error(key)
    iDelta == key && return property_not_implemented_error(key)
    iDmolar == key && return molar_density(model,result)
    iHmolar == key && return enthalpy(model,result)
    iSmolar == key && return entropy(model,result)
    iCpmolar == key && return isobaric_heat_capacity(model,result)
    iCp0molar == key && return isobaric_heat_capacity(idealmodel(model),result)
    iCvmolar == key && return isochoric_heat_capacity(model,result)
    iUmolar == key && return internal_energy(model,result)
    iGmolar == key && return gibbs_free_energy(model,result)
    iHelmholtzmolar == key && return helmholtz_energy(model,result)
    iHmolar_residual == key && return enthalpy_res(model,result)
    iSmolar_residual == key && return entropy_res(model,result)
    iGmolar_residual == key && return gibbs_free_energy_res(model,result)
    iHmolar_idealgas == key && return enthalpy(Clapeyron.idealmodel(model),result)
    iSmolar_idealgas == key && return entropy(Clapeyron.idealmodel(model),result)
    iUmolar_idealgas == key && return internal_energy(Clapeyron.idealmodel(model),result)
    iDmass == key && return mass_density(model,result)
    iHmass == key && return mass_enthalpy(model,result)
    iSmass == key && return mass_entropy(model,result)
    iCpmass == key && return mass_isobaric_heat_capacity(model,result)
    iCp0mass == key && return mass_isobaric_heat_capacity(idealmodel(model),result)
    iCvmass == key && return mass_isochoric_heat_capacity(model,result)
    iUmass == key && return mass_internal_energy(model,result)
    iGmass == key && return mass_gibbs_energy(model,result)
    iHelmholtzmass == key && return mass_helmholtz_energy(model,result)
    iHmass_idealgas == key && return mass_enthalpy(Clapeyron.idealmodel(model),result)
    iSmass_idealgas == key && return mass_entropy(Clapeyron.idealmodel(model),result)
    iUmass_idealgas == key && return mass_internal_energy(Clapeyron.idealmodel(model),result)
    iviscosity == key && return property_not_implemented_error(key)
    iconductivity == key && return property_not_implemented_error(key)
    isurface_tension == key && return property_not_implemented_error(key)
    iPrandtl == key && return property_not_implemented_error(key)
    ispeed_sound == key && return speed_of_sound(model,result)
    iisothermal_compressibility == key && return isothermal_compressibility(model,result)
    iisobaric_expansion_coefficient == key && return isobaric_expansivity(model,result)
    iisentropic_expansion_coefficient == key && return isentropic_compressibility(model,result)
    ifundamental_derivative_of_gas_dynamics == key && fundamental_derivative_of_gas_dynamics(model,result)
    ialphar == key && return return property_not_implemented_error(key)
    idalphar_dtau_constdelta == key && return property_not_implemented_error(key)
    idalphar_ddelta_consttau == key && return property_not_implemented_error(key)
    ialpha0 == key && return property_not_implemented_error(key)
    idalpha0_dtau_constdelta == key && return property_not_implemented_error(key)
    idalpha0_ddelta_consttau == key && return property_not_implemented_error(key)
    id2alpha0_ddelta2_consttau == key && return property_not_implemented_error(key)
    id3alpha0_ddelta3_consttau == key && return property_not_implemented_error(key)
    iBvirial == key && return property_not_implemented_error(key)
    iCvirial == key && return property_not_implemented_error(key)
    idBvirial_dT == key && return property_not_implemented_error(key)
    idCvirial_dT == key && return property_not_implemented_error(key)
    iZ == key && return compressibility_factor(model,result)
    iPIP == key && return Clapeyron.pip(model,result)
    ifraction_min == key && return property_not_implemented_error(key)
    ifraction_max == key && return property_not_implemented_error(key)
    iT_freeze == key && return property_not_implemented_error(key)
    iGWP20 == key && return property_not_implemented_error(key)
    iGWP100 == key && return property_not_implemented_error(key)
    iGWP500 == key && return property_not_implemented_error(key)
    iFH == key && return property_not_implemented_error(key)
    iHH == key && return property_not_implemented_error(key)
    iPH == key && return property_not_implemented_error(key)
    iODP == key && return property_not_implemented_error(key)
    iPhase == key && return property_not_implemented_error(key)
    iundefined_parameter == key &&  return property_not_implemented_error(key)
    end
end