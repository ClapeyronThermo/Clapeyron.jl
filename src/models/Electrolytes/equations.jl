function ElectrolyteFractionVector(model::ElectrolyteModel,M,z=[1.])
    ν = model.stoic_coeff
    Mw = model.puremodel.params.Mw.values.*1e-3
    isalts = 1:length(model.stoic_coeff)
    iions = 1:length(model.stoic_coeff[1])
    x_solv = z ./ (1+sum(z[j]*Mw[j] for j in model.isolvents)*(sum(M[k]*sum(ν[k][i] for i ∈ iions) for k ∈ isalts)))
    x_ions = [sum(M[k]*ν[k][l] for k ∈ isalts) / (1/sum(z[j]*Mw[j] for j in model.isolvents)+(sum(M[k]*sum(ν[k][i] for i ∈ iions) for k ∈ isalts))) for l ∈ iions]
    return append!(x_solv,x_ions)
end

function x0_volume(model::ElectrolyteModel,p,T,z; phase = :unknown)
    phase = Symbol(phase)
    if phase === :unknown || is_liquid(phase)
        return 1.5*x0_volume_liquid(model.puremodel,T,z)
    elseif is_vapour(phase)
        return 1.5*x0_volume_gas(model.puremodel,p,T,z)
    elseif is_supercritical(phase)
     else
        error("unreachable state on x0_volume")
    end
end