#pressure(model::ElectrolyteModel,V,T) = pressure(model,V,T,[0.5,0.5])

function ElectrolyteFractionVector(model::ElectrolyteModel,M,z=[1.])
    ν = model.stoic_coeff
    Mw = model.puremodel.params.Mw.values.*1e-3
    isalts = 1:length(model.stoic_coeff)
    iions = 1:length(model.stoic_coeff[1])
    x_solv = z ./ (1+sum(z[j]*Mw[j] for j in model.isolvents)*(sum(M[k]*sum(ν[k][i] for i ∈ iions) for k ∈ isalts)))
    x_ions = [sum(M[k]*ν[k][l] for k ∈ isalts) / (1/sum(z[j]*Mw[j] for j in model.isolvents)+(sum(M[k]*sum(ν[k][i] for i ∈ iions) for k ∈ isalts))) for l ∈ iions]
    return append!(x_solv,x_ions)
end

function FractionSalt(model::ElectrolyteModel,z)
    salts = model.salts
    isolv = model.isolvents
    isalts = model.isalts
    ν = model.stoic_coeff[:,1:length(salts)]

    z_new = z[isolv]
    z_ion  = z[isalts]

    z_salts = (z_ion'*inv(ν))'

    append!(z_new,z_salts)
    return z_new
end

function salt_stoichiometry(model::ElectrolyteModel,salts)
    ions = model.components[model.charge.!=0]
    ν = zeros(length(salts),length(ions))
    for i in 1:length(salts)
        v = salts[i][2]
        for j in 1:length(v)
            ν[i,v[j][1].==ions] .= v[j][2]
        end
        if sum(ν[i,:].*model.charge[model.charge.!=0])!==0.
            throw(ArgumentError("The salt $i is not electroneutral"))
        end
    end
    if any(sum(ν;dims=1).==0)
        throw(ArgumentError("Not all ions are involved in the salts"))
    end

    return ν
end

function molality_to_composition(model::ElectrolyteModel,salts,m,zsolv=[1.])
    ν = salt_stoichiometry(model,salts)

    ions = model.components[model.charge.!=0]
    neutral = model.components[model.charge.==0]
    Mw = mw(model.neutralmodel).*1e-3

    isalts = 1:length(salts)
    iions = 1:length(ions)
    ineutral = 1:length(neutral)
    x_solv = zsolv ./ (1+sum(zsolv[j]*Mw[j] for j in ineutral)*(sum(m[k]*sum(ν[k,i] for i ∈ iions) for k ∈ isalts)))
    x_ions = [sum(m[k]*ν[k,l] for k ∈ isalts) / (1/sum(zsolv[j]*Mw[j] for j in ineutral)+(sum(m[k]*sum(ν[k,i] for i ∈ iions) for k ∈ isalts))) for l ∈ iions]
    return append!(x_solv,x_ions)
end

export molality_to_composition


#=
function Obj_Sat(model::ElectrolyteModel, F, T, V_l, V_v,scales)
    ν = model.stoic_coeff[1]
    fun(_V) = eos(model, _V, T,[0.5,0.5])
    A_l,Av_l = Solvers.f∂f(fun,V_l)
    A_v,Av_v =Solvers.f∂f(fun,V_v)
    μ_l = VT_chemical_potential(model,V_l,T,[0.5,0.5])
    μ_v = VT_chemical_potential(model,V_v,T,[0.5,0.5])
    (p_scale,μ_scale) = scales
    F[1] = -(Av_l-Av_v)*p_scale
    F[2] = sum((μ_l-μ_v).*ν)*μ_scale
    return F
end =#

#=
function check_valid_sat_pure(model::ElectrolyteModel,P_sat,V_l,V_v,T,ε0 = 5e7)
    ε = abs(V_l-V_v)/(eps(typeof(V_l-V_v)))
    ε <= ε0 && return false
    _,dpdvl = p∂p∂V(model,V_l,T,[0.5,0.5])
    _,dpdvv = p∂p∂V(model,V_v,T,[0.5,0.5])
    return (dpdvl <= 0) && (dpdvv <= 0)
    #if ΔV > ε then Vl and Vv are different values
end =#

