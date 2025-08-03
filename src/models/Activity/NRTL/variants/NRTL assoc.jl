struct NRTLAssocParam <: EoSParam
    a::PairParam{Float64}
    b::PairParam{Float64}
    c::PairParam{Float64}
    Mw::SingleParam{Float64}
    δA::SingleParam{Float64} # ToDo: Add option to have multiple types of association site on the same component
    δD::SingleParam{Float64}
    nA::SingleParam{Float64}
    nD::SingleParam{Float64}
    rI::SingleParam{Float64}
end

abstract type NRTLAssocModel <: ActivityModel end

struct NRTLAssoc{c<:EoSModel} <: NRTLAssocModel
    components::Array{String,1}
    params::NRTLAssocParam
    puremodel::EoSVectorParam{c}
    references::Array{String,1}
end

export NRTLAssoc

NRTLAssoc

default_locations(::Type{NRTLAssoc}) = [
    "properties/molarmass.csv",
    "Activity/NRTL/NRTL_assoc_unlike.csv",
    "Activity/NRTL/NRTL_assoc.csv"
]

function NRTLAssoc(components; puremodel=PR,
    userlocations = String[], 
    pure_userlocations = String[],
    verbose = false,
    reference_state = nothing)

    formatted_components = format_components(components)
    params = getparams(
        formatted_components,
        default_locations(NRTLAssoc);
        userlocations = userlocations,
        asymmetricparams=["a","b","tau","alpha"],
        ignore_missing_singleparams=["a","b","Mw","tau","alpha"],
        verbose = verbose
    )

    if !__ismissing(params,"tau") && __ismissing(params,"a") && __ismissing(params,"b")
        a = params["tau"]
        b = PairParam("b",formatted_components)
    elseif __ismissing(params,"tau") 
        a = params["a"]
        b = params["b"]
    else
        throw(ArgumentError("NRTL: tau and (a,b) are mutually exclusive parameters, please provide only one of them"))
    end

    if !__ismissing(params,"alpha") && __ismissing(params,"c")
        c = params["alpha"]
    elseif __ismissing(params,"alpha")
        c = params["c"]
    else
        throw(ArgumentError("NRTL: `alpha` and `c` are mutually exclusive parameters, please provide only one of them"))
    end

    Mw  = get(params,"Mw",SingleParam("Mw",formatted_components))

    δA  = get(params, "δA", SingleParam("δA", formatted_components, fill(0.0, length(formatted_components))))
    δD  = get(params, "δD", SingleParam("δd", formatted_components, fill(0.0, length(formatted_components))))
    nA  = get(params, "nA", SingleParam("nA", formatted_components, fill(0.0, length(formatted_components))))
    nD  = get(params, "nD", SingleParam("nD", formatted_components, fill(0.0, length(formatted_components))))
    rI  = get(params, "rI", SingleParam("rI", formatted_components, fill(0.0, length(formatted_components))))

    _puremodel = init_puremodel(puremodel,components,pure_userlocations,verbose)
    packagedparams = NRTLAssocParam(a,b,c,Mw,δA,δD,nA,nD,rI)
    #                    Main NRTL article       Main Association NRTL article 
    references = String["10.1002/aic.690140124", "10.1002/aic.17061"]
    model = NRTLAssoc(formatted_components,packagedparams,_puremodel,references)
    set_reference_state!(model,reference_state,verbose = verbose)
    return model
end

function excess_g_assoc(model::NRTLAssocModel, p, T, z)
    Tz = Base.promote_eltype(T, z)
    R = R̄
    x = z ./ sum(z)
    comps = @comps
    n = sum(z)

    δA = model.params.δA.values
    δD = model.params.δD.values
    nA = model.params.nA.values
    nD = model.params.nD.values
    rI = model.params.rI.values 

    ρA_i = nA ./ rI
    ρD_i = nD ./ rI

    ρA = sum(x .* ρA_i)
    ρD = sum(x .* ρD_i)

    Κref = Tz(0.034)
    εref = Tz(1960.0)
    Δref = Κref * (exp(εref / T) - one(Tz))

    ΔAD = δA .* δD' .* Δref

    XA = ones(Tz, length(comps))
    XD = ones(Tz, length(comps))

    for _ in 1:100
        XA_old = copy(XA)
        XD_old = copy(XD)
        for i ∈ 1:length(comps)
            denom_A = one(Tz) + sum(x[j] * ρD_i[j] * ΔAD[i,j] * XD[j] for j ∈ 1:length(comps))
            denom_D = one(Tz) + sum(x[j] * ρA_i[j] * ΔAD[j,i] * XA[j] for j ∈ 1:length(comps))
            XA[i] = one(Tz) / denom_A
            XD[i] = one(Tz) / denom_D
        end
        if maximum(abs.(XA .- XA_old)) < 1e-10 && maximum(abs.(XD .- XD_old)) < 1e-10
            break
        end
    end

    lnγ_assoc = zero(Tz)
    for i ∈ 1:length(comps)
        termA = nA[i] * (log(XA[i]) + Tz(0.5) * (XA[i] - one(Tz)))
        termD = nD[i] * (log(XD[i]) + Tz(0.5) * (XD[i] - one(Tz)))
        lnγ_assoc += x[i] * (termA + termD)
    end

    return n * R * T * lnγ_assoc
end

function excess_g_comb(model::NRTLAssocModel,p,T,z)
    Tz = Base.promote_eltype(T, z)
    R = R̄
    x = z ./ sum(z)
    comps = @comps
    n = sum(z)

    rI = model.params.rI.values

    ϕ_I = rI ^ (2/3) .* x ./ sum(rI ^ (2/3) .* x)
    lnγ_comb = one(Tz) - ϕ_I./x - ln(ϕ_I./x)

    return n * R * T * sum(lnγ_comb)
end

function excess_g_res(model::NRTLAssocModel,p,T,z)
    a = model.params.a.values
    b  = model.params.b.values
    c  = model.params.c.values
    _0 = zero(Base.promote_eltype(model,p,T,z))
    n = sum(z)
    invn = 1/n
    invT = 1/(T)
    res = _0 
    for i ∈ @comps
        ΣτGx = _0
        ΣGx = _0
        xi = z[i]*invn
        for j ∈ @comps
            xj = z[j]*invn
            τji = a[j,i] + b[j,i]*invT
            Gji = exp(-c[j,i]*τji)
            Gx = xj*Gji
            ΣGx += Gx
            ΣτGx += Gx*τji
        end
        res += xi*ΣτGx/ΣGx
    end
    return n*res*R̄*T

end

function excess_gibbs_free_energy(model::NRTLAssocModel,p,T,z) 

    g_assoc = excess_g_assoc(model, p, T, z)
    g_comb = excess_g_comb(model, p, T, z)
    g_res = excess_g_res(model,p,T,z)
    return g_comb + g_res + g_assoc 

end    