struct COSMOSAC02Param <: EoSParam
    Pi::SingleParam{Vector{Float64}}
    V::SingleParam{Float64}
    A::SingleParam{Float64}
end

abstract type COSMOSAC02Model <: ActivityModel end

struct COSMOSAC02{c<:EoSModel} <: COSMOSAC02Model
    components::Array{String,1}
    params::COSMOSAC02Param
    puremodel::EoSVectorParam{c}
    absolutetolerance::Float64
    references::Array{String,1}
end

export COSMOSAC02
"""
    COSMOSAC02(components;
    puremodel = PR,
    userlocations = String[],
    pure_userlocations = String[],
    verbose = false)

## Input parameters:
- `Pi` :Single Parameter{String} 
- `V`: Single Parameter{Float64}
- `A`: Single Parameter{Float64}

## Description
An activity coefficient model using molecular solvation based on the COSMO-RS method.

## References
1. Lin, S-T. & Sandler, S.I. (2002). A priori phase equilibrium prediction from a segment contribution solvation model. Industrial & Engineering Chemistry Research, 41(5), 899–913. [doi:10.1021/ie001047w](https://doi.org/10.1021/ie001047w)
"""
COSMOSAC02

default_locations(::Type{COSMOSAC02}) = ["Activity/COSMOSAC/COSMOSAC02_like.csv"]

function COSMOSAC02(components;
    puremodel = PR,
    userlocations = String[],
    pure_userlocations = String[],
    use_nist_database = false,
    verbose = false)
    formatted_components = format_components(components)

    if use_nist_database
        @warn "using parameters from the nistgov/COSMOSAC database, check their license before usage."
        CAS, INCHIKEY = get_cosmo_comps()
        A = zeros(length(components))
        V = zeros(length(components))
        Pi = [zeros(51) for i in 1:length(components)]
        for i in 1:length(components)
            id = cas(formatted_components[i])
            ids = CAS.==uppercase(id[1])
            dbname = INCHIKEY[ids]
            file = String(take!(Downloads.download("https://raw.githubusercontent.com/usnistgov/COSMOSAC/master/profiles/UD/sigma/"*dbname[1]*".sigma", IOBuffer())))
            lines = split(file,r"\n")
            meta = lines[1][9:end]
            json = JSON3.read(meta)
            A[i] = json["area [A^2]"]
            V[i] = json["volume [A^3]"]
            Pi[i] = [parse(Float64,split(lines[i]," ")[2]) for i in 4:54]
        end
        A = SingleParam("A",formatted_components,A)
        V = SingleParam("V",formatted_components,V)
        Pi = SingleParam("Pi",formatted_components,Pi)
    else
        params = getparams(formatted_components, default_locations(COSMOSAC02); userlocations = userlocations, verbose = verbose)
        Pi  = COSMO_parse_Pi(params["Pi"])
        A  = params["A"]
        V  = params["V"]
    end


    _puremodel = init_puremodel(puremodel,components,pure_userlocations,verbose)
    packagedparams = COSMOSAC02Param(Pi,V,A)
    references = String["10.1021/ie001047w","10.1021/acs.jctc.9b01016","10.1021/acs.iecr.7b01360"]
    model = COSMOSAC02(formatted_components,packagedparams,_puremodel,1e-12,references)
    return model
end

function activity_coefficient(model::COSMOSAC02Model,V,T,z)
    return exp.(@f(lnγ_comb) .+ @f(lnγ_res))
end

function excess_g_res(model::COSMOSAC02Model,V,T,z)
    lnγ = @f(lnγ_res)
    sum(z[i]*R̄*T*lnγ[i] for i ∈ @comps)
end

function lnγ_comb(model::COSMOSAC02Model,p,T,z)
    r0 = 66.69
    q0 = 79.53

    V = model.params.V.values
    A = model.params.A.values
    ∑z = sum(z)
    r̄⁻¹ = ∑z/dot(V,z)
    q̄⁻¹ = ∑z/dot(A,z)
    lnγ_comb = zeros(typeof(T+first(z)),length(z))
    for i ∈ @comps
        qi = A[i]/q0
        Φi = V[i]*r̄⁻¹
        θi = q0*qi*q̄⁻¹
        lnγ_comb[i] = log(Φi)+(1-Φi)-5*qi*(log(Φi/θi)+(1-Φi/θi))
    end
    return lnγ_comb
  end

function lnγ_res(model::COSMOSAC02Model,V,T,z)
    aeff = 7.5
    A = model.params.A.values
    n = A ./ aeff
    Pi = model.params.Pi.values
    #PS = sum(x[i]*Pi[i][:] for i ∈ @comps) ./ sum(x[i]*A[i] for i ∈ @comps)
    PS = zeros(eltype(z),length(Pi[1]))
    for i ∈ @comps
        PS .+= Ref(z[i]) .* Pi[i]       
    end
    PS .*= 1/dot(z,A)
    
    lnΓS = @f(lnΓ,PS)
    lnΓi = [@f(lnΓ,Pi[i]./A[i]) for i ∈ @comps]
    lnγ_res_ = [n[i]*sum(Pi[i][v]/A[i]*(lnΓS[v]-lnΓi[i][v]) for v ∈ 1:51) for i ∈ @comps]
    
    return lnγ_res_
end

function lnΓ(model::COSMOSAC02Model,V,T,z,P)
    Γold =lnΓ_fixpoint(P,T,z)
    atol = model.absolutetolerance
    lnΓ_f0(x,y) = lnΓ_fixpoint(x,y,P,T)
    Γ = Solvers.fixpoint(lnΓ_f0,Γold,Solvers.SSFixPoint(0.5),atol =atol ,max_iters=1000)
    Γ .= log.(Γ)
    return Γ
end

function lnΓ_fixpoint(Γnew,Γold,P,T)
    Γnew .= zero(eltype(Γnew))
    σ  = -0.025:0.001:0.025
    Tinv = one(T)/T
    @inbounds for i = 1:length(Γnew)
        for j = 1:length(Γnew)
            Γnew[j] += P[i]*Γold[i]*exp(-ΔW(σ[j],σ[i])*Tinv) 
        end
    end
    Γnew .= one(eltype(Γnew))./Γnew
end

function lnΓ_fixpoint(P,T,z)
    Γnew = zeros(promote_type(eltype(T),eltype(P),eltype(z)),length(P))
    σ  = -0.025:0.001:0.025
    Tinv = one(T)/T
    @inbounds for i = 1:length(Γnew)
        for j = 1:length(Γnew)
            Γnew[j] += P[i]*exp(-ΔW(σ[j],σ[i])*Tinv) 
        end
    end
    Γnew .= one(eltype(Γnew))./Γnew
end

function ΔW(σm,σn)
    σdon,σacc = minmax(σm,σn)
    σhb  = 0.0084
    chb  = 85580
    α    = 16466.72
    R    = 0.001987
    return (α/2*(σm+σn)^2+chb*max(0,σacc-σhb)*min(0,σdon+σhb))/R
end


#=
function lnγ_comb(model::COSMOSAC02Model,V,T,z)
    r0 = 66.69
    q0 = 79.53

    V = model.params.V.values
    A = model.params.A.values

    x = z ./ sum(z)

    r = V ./ r0
    q = A ./ q0

    Φ = r/sum(x[i]*r[i] for i ∈ @comps)
    θ = q/sum(x[i]*q[i] for i ∈ @comps)
    lnγ_comb = @. log(Φ)+(1-Φ)-5*q*(log(Φ/θ)+(1-Φ/θ))
    return lnγ_comb
end
=#
