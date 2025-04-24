function QT_property(model,q,T,z,f::F,p0) where F
    if f == temperature
        return T
    end

    res = qt_flash(model,q,T,z,p0 = p0)
    if f == temperature
        return temperature(res)
    elseif p == pressure
        return pressure(res)
    else
        return f(model,res)
    end
end

module QT
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
        function $f(model,q,T,z = Clapeyron.SA[1.0],p0 = nothing)
            Clapeyron.QT_property(model,q,T,z,Clapeyron.$f,p0)
        end
    end

    function flash(model,q,T,z = Clapeyron.SA[1.0],args...;kwargs...)
        return Clapeyron.qt_flash(model,q,T,z,args...;kwargs...)
    end
end
#export chemical_potential, activity_coefficient, activity, aqueous_activity, fugacity_coefficient,reference_chemical_potential,reference_chemical_potential_type
#export chemical_potential_res
#export mixing, excess, gibbs_solvation

end #module

function QP_property(model,q,p,z,f::F,T0) where F
    if f == pressure
        return p
    end

    res = qp_flash(model,q,p,z,T0 = T0)
    if f == temperature
        return temperature(res)
    elseif p == pressure
        return pressure(res)
    else
        return f(model,res)
    end
end

module QP
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
 :mass_density, :molar_density, :compressibility_factor,
#molar :gradient :properties
 :identify_phase]
    @eval begin
        function $f(model,q,p,z = Clapeyron.SA[1.0],T0 = nothing)
            Clapeyron.QP_property(model,q,p,z,Clapeyron.$f,T0)
        end
    end

    function flash(model,q,p,z = Clapeyron.SA[1.0],args...;kwargs...)
        return Clapeyron.qp_flash(model,q,p,z,args...;kwargs...)
    end
end
#export chemical_potential, activity_coefficient, activity, aqueous_activity, fugacity_coefficient,reference_chemical_potential,reference_chemical_potential_type
#export chemical_potential_res
#export mixing, excess, gibbs_solvation

end #module