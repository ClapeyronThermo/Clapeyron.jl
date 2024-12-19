function PH_property(model,p,h,z,f::F,phase,T0,threaded) where F
    if f == enthalpy
        return h
    end

    if f == pressure
        return p
    end

    if !is_unknown(phase)
        T,calc_phase = _Tproperty(model,p,h,z,T0 = T0,phase = phase,threaded = threaded)
        if calc_phase != :eq && calc_phase != :failure
            return f(model,p,T,z;phase = calc_phase)
        elseif calc_phase == :eq && !supports_lever_rule(f)
            thow(invalid_property_multiphase_error(f))
        elseif calc_phase == :eq && supports_lever_rule(f)
            result = ph_flash(model,p,h,z,T0)
            return f(model,result)
        else
            return f(model,p,T,z;phase = phase)
        end
    else
        res = ph_flash(model,p,h,z,T0 = T0)
        if f == temperature
            return temperature(res)
        else
            return f(model,res)
        end
    end
end

module PH
import Clapeyron
for f in [:temperature,:volume, :pressure, :entropy, :internal_energy, :enthalpy, :gibbs_free_energy, :helmholtz_free_energy,
 :entropy_res, :internal_energy_res, :enthalpy_res, :gibbs_free_energy_res, :helmholtz_free_energy_res,
#second derivative order properties
 :isochoric_heat_capacity, :isobaric_heat_capacity, :adiabatic_index,
 :isothermal_compressibility, :isentropic_compressibility, :speed_of_sound,
 :isobaric_expansivity, :joule_thomson_coefficient, :inversion_temperature,
#higher :derivative :order :properties
 :fundamental_derivative_of_gas_dynamics,
#volume :properties
 :mass_density,:molar_density, :compressibility_factor,
#molar :gradient :properties
 :identify_phase]
    @eval begin
        function $f(model,p,h,z = Clapeyron.SA[1.0];phase = :unknown,T0 = nothing, threaded = true)
            Clapeyron.PH_property(model,p,h,z,Clapeyron.$f,phase,T0,threaded)
        end
    end
end
#export chemical_potential, activity_coefficient, activity, aqueous_activity, fugacity_coefficient,reference_chemical_potential,reference_chemical_potential_type
#export chemical_potential_res
#export mixing, excess, gibbs_solvation

end  #module