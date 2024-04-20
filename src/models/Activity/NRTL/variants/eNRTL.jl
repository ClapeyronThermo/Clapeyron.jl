abstract type eNRTLModel <: ElectrolyteActivityModel end

struct eNRTL{T<:IdealModel,c<:EoSModel,i<:IonModel} <: eNRTLModel
    components::Array{String,1}
    icomponents::UnitRange{Int}
    charge::Vector{Int64}
    idealmodel::T
    neutralmodel::c
    ionmodel::i
    references::Array{String,1}
end
function eNRTL(solvents,ions; 
    idealmodel = BasicIdeal,
    neutralmodel = NRTL,
    ionmodel = PDH,
    RSPmodel = ConstRSP,
    userlocations=String[], 
    ideal_userlocations=String[],
     verbose=false)
    components = deepcopy(ions)
    prepend!(components,solvents)

    params = getparams(components, ["Electrolytes/properties/charges.csv"]; userlocations=userlocations, verbose=verbose)
    charge = params["charge"].values

    icomponents = 1:length(components)

    neutral_path = [DB_PATH*"/"*default_locations(neutralmodel)[1],DB_PATH*"/Activity/NRTL/eNRTL"]

    init_idealmodel = init_model(idealmodel,components,ideal_userlocations,verbose)
    init_neutralmodel = neutralmodel(components;userlocations=append!(userlocations,neutral_path),verbose=verbose)
    init_ionmodel = ionmodel(solvents,ions;RSPmodel=RSPmodel,userlocations=append!(userlocations,neutral_path),verbose=verbose)

    references = String[]
    model = eNRTL(components,icomponents,charge,init_idealmodel,init_neutralmodel,init_ionmodel,references)
    return model
end

export eNRTL

function excess_g_res(model::eNRTLModel, V, T, z)

    return excess_g_local(model,V,T,z)+excess_g_res(model.ionmodel,V,T,z)
end

function excess_g_local(model::eNRTLModel, V, T, z)
    τ = model.neutralmodel.params.a.values
    b  = model.neutralmodel.params.b.values
    α  = model.neutralmodel.params.c.values
    _0 = zero(T+first(z))
    n = sum(z)
    invn = 1/n
    invT = 1/(T)
    res = _0 


    isolv = model.icomponents[model.charge.==0]
    icat = model.icomponents[model.charge.>0]
    iani = model.icomponents[model.charge.<0]
    
    Ya = z[iani]
    Ya ./= sum(Ya)

    Yc = z[icat]
    Yc ./= sum(Yc)

    for m in isolv
        ∑τGx = _0
        ∑Gx = _0
        xm = z[m]*invn
        for j in @comps
            xj = z[j]*invn
            τjm = τ[j,m] + b[j,m]*invT
            Gjm = exp(-α[j,m]*τjm)
            Gx = xj*Gjm
            ∑Gx += Gx
            ∑τGx += Gx*τjm
        end
        res += xm*∑τGx/∑Gx
    end

    for c in icat
        ∑τGx = _0
        ∑Gx = _0
        xc = z[c]*invn
        Zc = abs(model.charge[c])
        for j in @comps
            xj = z[j]*invn
            τjc = τ[j,c] + b[j,c]*invT
            Gjc = exp(-α[j,c]*τjc)
            Gx = xj*Gjc
            ∑Gx += Gx
            ∑τGx += Gx*τjc
        end
        res += xc*Zc*∑τGx/∑Gx
    end

    for a in iani
        ∑τGx = _0
        ∑Gx = _0
        xa = z[a]*invn
        Za = abs(model.charge[a])
        for j in @comps
            xj = z[j]*invn
            τja = τ[j,a] + b[j,a]*invT
            Gja = exp(-α[j,a]*τja)
            Gx = xj*Gja
            ∑Gx += Gx
            ∑τGx += Gx*τja
        end
        res += xa*Za*∑τGx/∑Gx
    end
    return n*res*R̄*T
end

function activity_coefficient(model::eNRTLModel,p,T,z)
    ion = model.charge.!=0

    X = gradient_type(p,T,z)
    lnγres = (Solvers.gradient(x->excess_g_local(model,p,T,x),z)/(R̄*T))::X
    lnγion = (Solvers.gradient(x->excess_gibbs_free_energy(model.ionmodel,p,T,x),z)/(R̄*T))::X

    z0 = ones(length(z))*1e-30
    @. z0[!ion] = z[!ion]
    z0 = z0./sum(z0)

    lnγres0 = (Solvers.gradient(x->excess_g_local(model,p,T,x),z0)/(R̄*T))::X


    lnγ = @. ion*(lnγres-lnγres0)+!ion*lnγres+lnγion

    return exp.(lnγ)
end