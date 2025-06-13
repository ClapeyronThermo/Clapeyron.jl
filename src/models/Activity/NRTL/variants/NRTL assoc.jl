struct NRTLAssocParam <: EoSParam
    a::PairParam{Float64}
    b::PairParam{Float64}
    c::PairParam{Float64}
    Mw::SingleParam{Float64}
    δA::SingleParam{Float64}
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

default_locations(::Type{NRTLAssoc}) = ["properties/molarmass.csv","Activity/NRTL/NRTL_unlike.csv","Activity/NRTL/NRTL_assoc.csv"]

function NRTLAssoc(components; puremodel=PR,
    userlocations = String[], 
    pure_userlocations = String[],
    verbose = false,
    reference_state = nothing)

    formatted_components = format_components(components)
    params = getparams(formatted_components, default_locations(NRTLAssoc); userlocations = userlocations, asymmetricparams=["a","b","tau","alpha"], ignore_missing_singleparams=["a","b","Mw","tau","alpha"], verbose = verbose)
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
    # Association parameters with safe defaults
    δA  = haskey(params, "δA") ? params["δA"] : SingleParam("δA", formatted_components; values = fill(0.0, length(formatted_components)))
    δD  = haskey(params, "δD") ? params["δD"] : SingleParam("δD", formatted_components; values = fill(0.0, length(formatted_components)))
    nA  = haskey(params, "nA") ? params["nA"] : SingleParam("nA", formatted_components; values = fill(0.0, length(formatted_components)))
    nD  = haskey(params, "nD") ? params["nD"] : SingleParam("nD", formatted_components; values = fill(0.0, length(formatted_components)))
    rI  = haskey(params, "rI") ? params["rI"] : SingleParam("rI", formatted_components; values = fill(1.0, length(formatted_components)))
    
    _puremodel = init_puremodel(puremodel,components,pure_userlocations,verbose)
    packagedparams = NRTLAssocParam(a,b,c,Mw,δA,δD,nA,nD,rI)
    references = String["10.1002/aic.690140124"]
    model = NRTLAssoc(formatted_components,packagedparams,_puremodel,references)
    set_reference_state!(model,reference_state,verbose = verbose)
    return model
end

#=
function activity_coefficient(model::NRTLModel,p,T,z)
    a = model.params.a.values
    b = model.params.b.values
    c = model.params.c.values
    x = z ./ sum(z)
    τ = @. a+b/T
    G = @. exp(-c*τ)
    lnγ = sum(x[j]*τ[j,:].*G[j,:] for j ∈ @comps)./sum(x[k]*G[k,:] for k ∈ @comps)+sum(x[j]*G[:,j]/sum(x[k]*G[k,j] for k ∈ @comps).*(τ[:,j] .-sum(x[m]*τ[m,j]*G[m,j] for m ∈ @comps)/sum(x[k]*G[k,j] for k ∈ @comps)) for j in @comps)
    return exp.(lnγ)
end
=#

function excess_g_assoc(model::NRTLAssocModel, p, T, z)
    R = R̄
    x = z ./ sum(z)
    comps = @comps
    n = sum(z)

    δA = model.params.δA.values
    δD = model.params.δD.values
    nA = model.params.nA.values
    nD = model.params.nD.values
    rI = model.params.rI.values 

    # Define site types
    site_types = [:A, :D]

    # ρA_i = nA_i / r_i
    ρA_i = nA ./ rI
    ρD_i = nD ./ rI

    # ρA = ∑ xᵢ * nAᵢ / rᵢ
    ρA = sum(x .* ρA_i)
    ρD = sum(x .* ρD_i)

    # Compute reference Δ_AD 
    κref = 0.034
    εref = 1960.0
    Δref = κref * (exp(εref / T) - 1)

    # Construct Δ_AB matrix (assuming only A-D interaction matters)
    ΔAD = δA .* δD' .* Δref

    # Site fractions: solve fixed-point for X_A and X_D
    XA = one(eltype(z)) .* ones(length(comps))
    XD = one(eltype(z)) .* ones(length(comps))

    # Iterative solve for X_A
    for _ in 1:100
        XA_old = copy(XA)
        XD_old = copy(XD)
        for i ∈ 1:length(comps)
            denom_A = 1.0 + sum(x[j] * ρD_i[j] * ΔAD[i,j] * XD[j] for j ∈ 1:length(comps))
            denom_D = 1.0 + sum(x[j] * ρA_i[j] * ΔAD[j,i] * XA[j] for j ∈ 1:length(comps))
            XA[i] = 1.0 / denom_A
            XD[i] = 1.0 / denom_D
        end
        if maximum(abs.(XA .- XA_old)) < 1e-10 &&
           maximum(abs.(XD .- XD_old)) < 1e-10
            break
        end
    end

    # Now compute gE_assoc 
    lnγ_assoc = zero(eltype(z))
    for i ∈ 1:length(comps)
        termA = nA[i] * (log(XA[i]) + 0.5 * (XA[i] - 1.0))
        termD = nD[i] * (log(XD[i]) + 0.5 * (XD[i] - 1.0))
        lnγ_assoc += x[i] * (termA + termD)
    end

    return n * R * T * lnγ_assoc
end

function excess_g_res(model::NRTLAssocModel,p,T,z)
    g_assoc = excess_g_assoc(model, p, T, z)
    a = model.params.a.values
    b  = model.params.b.values
    c  = model.params.c.values
    _0 = zero(Base.promote_eltype(model,p,T,z))
    n = sum(z)
    invn = 1/n
    invT = 1/(T)
    res = _0 
    for i ∈ @comps
        ∑τGx = _0
        ∑Gx = _0
        xi = z[i]*invn
        for j ∈ @comps
            xj = z[j]*invn
            τji = a[j,i] + b[j,i]*invT
            Gji = exp(-c[j,i]*τji)
            Gx = xj*Gji
            ∑Gx += Gx
            ∑τGx += Gx*τji
        end
        res += xi*∑τGx/∑Gx
    end
    return n*res*R̄*T
end

excess_gibbs_free_energy(model::NRTLAssocModel,p,T,z) = excess_g_res(model,p,T,z)