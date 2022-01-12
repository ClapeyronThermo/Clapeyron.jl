struct COSMOSAC02Param <: EoSParam
    Pi::SingleParam{Vector{Float64}}
    V::SingleParam{Float64}
    A::SingleParam{Float64}
end

abstract type COSMOSAC02Model <: ActivityModel end

struct COSMOSAC02{c<:EoSModel} <: COSMOSAC02Model
    components::Array{String,1}
    icomponents::UnitRange{Int}
    params::COSMOSAC02Param
    puremodel::Vector{c}
    absolutetolerance::Float64
    references::Array{String,1}
end

@registermodel COSMOSAC02
export COSMOSAC02

function COSMOSAC02(components; puremodel=PR,
    userlocations=String[], 
     verbose=false)
    params = getparams(components, ["Activity/COSMOSAC/COSMOSAC02_like.csv"]; userlocations=userlocations, verbose=verbose)
    Pi  = COSMO_parse_Pi(params["Pi"])
    A  = params["A"]
    V  = params["V"]
    icomponents = 1:length(components) 


    init_puremodel = [puremodel([components[i]]) for i in icomponents]
    packagedparams = COSMOSAC02Param(Pi,V,A)
    references = String[]
    model = COSMOSAC02(components,icomponents,packagedparams,init_puremodel,1e-12,references)
    return model
end

function activity_coefficient(model::COSMOSAC02Model,V,T,z)
    return exp.(@f(lnγ_comb) .+ @f(lnγ_res))
end

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

function lnγ_res(model::COSMOSAC02Model,V,T,z)
    x = z ./ sum(z)

    aeff = 7.5
    A = model.params.A.values
    n = A ./ aeff
    Pi = model.params.Pi.values
    PS = sum(x[i]*Pi[i][:] for i ∈ @comps) ./ sum(x[i]*A[i] for i ∈ @comps)

    lnΓS = @f(lnΓ,PS)
    lnΓi = [@f(lnΓ,Pi[i]./A[i]) for i ∈ @comps]
    lnγ_res_ =  [n[i]*sum(Pi[i][v]/A[i]*(lnΓS[v]-lnΓi[i][v]) for v ∈ 1:51) for i ∈ @comps]
    
    return lnγ_res_
end

function lnΓ(model::COSMOSAC02Model,V,T,z,P)
    Γ0 = ones(length(P))
    σ  = -0.025:0.001:0.025
    Γold = exp.(-log.(sum(P[i]*Γ0[i]*exp.(-ΔW.(σ,σ[i])./T) for i ∈ 1:51)))
    Γnew = Γold
    tol = 1
    i = 1
    damp_factor = 0.5
    while tol>sqrt(model.absolutetolerance)
        Γnew = exp.(-log.(sum(P[i]*Γold[i]*exp.(-ΔW.(σ,σ[i])./T) for i ∈ 1:51)))
        Γnew .*= (1-damp_factor)
        Γnew .+= damp_factor .* Γold
        tol = cosmo_tol(Γnew,Γold)
        Γold .= Γnew
        i+=1
    end
    lnΓ = log.(Γold)
    return lnΓ
end

function ΔW(σm,σn)
    σdon,σacc = minmax(σm,σn)
    σhb  = 0.0084
    chb  = 85580
    α    = 16466.72
    R    = 0.001987
    return (α/2*(σm+σn)^2+chb*max(0,σacc-σhb)*min(0,σdon+σhb))/R
end