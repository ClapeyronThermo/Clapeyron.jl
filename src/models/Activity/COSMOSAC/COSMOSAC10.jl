struct COSMOSAC10Param <: EoSParam
    Pnhb::SingleParam{Vector{Float64}}
    POH::SingleParam{Vector{Float64}}
    POT::SingleParam{Vector{Float64}}
    V::SingleParam{Float64}
    A::SingleParam{Float64}
end

abstract type COSMOSAC10Model <: COSMOSAC02Model end

struct COSMOSAC10{c<:EoSModel} <: COSMOSAC10Model
    components::Array{String,1}
    params::COSMOSAC10Param
    puremodel::EoSVectorParam{c}
    absolutetolerance::Float64
    references::Array{String,1}
end

@registermodel COSMOSAC10
export COSMOSAC10

function COSMOSAC10(components::Vector{String};
    puremodel = PR,
    userlocations = String[],
    pure_userlocations = String[],
    verbose=false, kwargs...)

    params = getparams(components, ["Activity/COSMOSAC/COSMOSAC10_like.csv"]; userlocations=userlocations, verbose=verbose)
    Pnhb  = COSMO_parse_Pi(params["Pnhb"])
    POH  = COSMO_parse_Pi(params["POH"])
    POT  = COSMO_parse_Pi(params["POT"])
    A  = params["A"]
    V  = params["V"]
    
    _puremodel = init_puremodel(puremodel,components,pure_userlocations,verbose)    
    packagedparams = COSMOSAC10Param(Pnhb,POH,POT,V,A)
    references = String[]
    model = COSMOSAC10(components,packagedparams,_puremodel,1e-12,references)
    return model
end

function lnγ_res(model::COSMOSAC10Model,V,T,z)
    x = z ./ sum(z)

    aeff = 7.5
    A = model.params.A.values
    n = A ./ aeff

    Pnhb  = model.params.Pnhb.values
    PSnhb = sum(x[i]*Pnhb[i] for i ∈ @comps) ./ sum(x[i]*A[i] for i ∈ @comps)

    POH = model.params.POH.values
    PSOH = sum(x[i]*POH[i] for i ∈ @comps) ./ sum(x[i]*A[i] for i ∈ @comps)

    POT = model.params.POT.values
    PSOT = sum(x[i]*POT[i] for i ∈ @comps) ./ sum(x[i]*A[i] for i ∈ @comps)

    (lnΓSnhb, lnΓSOH, lnΓSOT)= @f(lnΓ,PSnhb,PSOH,PSOT)
    lnΓi = [@f(lnΓ,Pnhb[i]./A[i],POH[i]./A[i],POT[i]./A[i]) for i ∈ @comps]
    
    lnγ_res_ =  [n[i]*(sum(Pnhb[i][v]/A[i]*(lnΓSnhb[v]-lnΓi[i][1][v]) for v ∈ 1:51)
                      +sum(POH[i][v]/A[i]*(lnΓSOH[v]-lnΓi[i][2][v]) for v ∈ 1:51)
                      +sum(POT[i][v]/A[i]*(lnΓSOT[v]-lnΓi[i][3][v]) for v ∈ 1:51)) for i ∈ @comps]               
    return lnγ_res_
end

function lnΓ(model::COSMOSAC10Model,V,T,z,Pnhb,POH,POT)
    _TYPE = typeof(V+T+first(z))
    Γ0nhb = ones(_TYPE,length(Pnhb))
    Γ0OH = ones(_TYPE,length(POH))
    Γ0OT = ones(_TYPE,length(POT))

    σ   = -0.025:0.001:0.025
    Γnhb_old = exp.(-log.(sum(Pnhb[i]*Γ0nhb[i]*exp.(-ΔW.(σ,σ[i],1,1,T)./T) for i ∈ 1:51)
                         +sum(POH[i]*Γ0OH[i]*exp.(-ΔW.(σ,σ[i],2,1,T)./T) for i ∈ 1:51)
                         +sum(POT[i]*Γ0OT[i]*exp.(-ΔW.(σ,σ[i],3,1,T)./T) for i ∈ 1:51)))
    ΓOH_old = exp.(-log.(sum(Pnhb[i]*Γ0nhb[i]*exp.(-ΔW.(σ,σ[i],1,2,T)./T) for i ∈ 1:51)
                         +sum(POH[i]*Γ0OH[i]*exp.(-ΔW.(σ,σ[i],2,2,T)./T) for i ∈ 1:51)
                         +sum(POT[i]*Γ0OT[i]*exp.(-ΔW.(σ,σ[i],3,2,T)./T) for i ∈ 1:51)))
    ΓOT_old = exp.(-log.(sum(Pnhb[i]*Γ0nhb[i]*exp.(-ΔW.(σ,σ[i],1,3,T)./T) for i ∈ 1:51)
                         +sum(POH[i]*Γ0OH[i]*exp.(-ΔW.(σ,σ[i],2,3,T)./T) for i ∈ 1:51)
                         +sum(POT[i]*Γ0OT[i]*exp.(-ΔW.(σ,σ[i],3,3,T)./T) for i ∈ 1:51)))
    Γnhb_new = deepcopy(Γnhb_old)
    ΓOH_new = deepcopy(ΓOH_old)
    ΓOT_new = deepcopy(ΓOT_old)
    tol = 1
    i = 1

    while tol>sqrt(model.absolutetolerance)
        Γnhb_new .= exp.(-log.(sum(Pnhb[i]*Γnhb_old[i]*exp.(-ΔW.(σ,σ[i],1,1,T)./T) for i ∈ 1:51)
                         +sum(POH[i]*ΓOH_old[i]*exp.(-ΔW.(σ,σ[i],2,1,T)./T) for i ∈ 1:51)
                         +sum(POT[i]*ΓOT_old[i]*exp.(-ΔW.(σ,σ[i],3,1,T)./T) for i ∈ 1:51)))
        ΓOH_new .= exp.(-log.(sum(Pnhb[i]*Γnhb_old[i]*exp.(-ΔW.(σ,σ[i],1,2,T)./T) for i ∈ 1:51)
                         +sum(POH[i]*ΓOH_old[i]*exp.(-ΔW.(σ,σ[i],2,2,T)./T) for i ∈ 1:51)
                         +sum(POT[i]*ΓOT_old[i]*exp.(-ΔW.(σ,σ[i],3,2,T)./T) for i ∈ 1:51)))
        ΓOT_new .= exp.(-log.(sum(Pnhb[i]*Γnhb_old[i]*exp.(-ΔW.(σ,σ[i],1,3,T)./T) for i ∈ 1:51)
                         +sum(POH[i]*ΓOH_old[i]*exp.(-ΔW.(σ,σ[i],2,3,T)./T) for i ∈ 1:51)
                         +sum(POT[i]*ΓOT_old[i]*exp.(-ΔW.(σ,σ[i],3,3,T)./T) for i ∈ 1:51)))

        Γnhb_new .= (Γnhb_new .+ Γnhb_old) ./2
        ΓOH_new .= (ΓOH_new .+ ΓOH_old) ./2
        ΓOT_new .= (ΓOT_new .+ ΓOT_old) ./2

        tol = (cosmo_tol(Γnhb_new,Γnhb_old) + cosmo_tol(ΓOH_new,ΓOH_old) + cosmo_tol(ΓOT_new,ΓOT_old))/3
        Γnhb_old .= Γnhb_new
        ΓOH_old .= ΓOH_new
        ΓOT_old .= ΓOT_new
        i+=1
    end
    lnΓnhb = log.(Γnhb_old)
    lnΓOH = log.(ΓOH_old)
    lnΓOT = log.(ΓOT_old)

    return (lnΓnhb,lnΓOH,lnΓOT)
end

function ΔW(σm,σn,t,s,T)
    ces  = 6525.69+1.4859e8/T^2
    if t==2 && s==2 && σm*σn<0
        chb = 4013.78
    elseif t==3 && s==3 && σm*σn<0
        chb = 932.31
    elseif t==3 && s==2 && σm*σn<0
        chb = 3016.43
    else
        chb = 0.
    end
    R    = 0.001987
    return (ces*(σm+σn)^2-chb*(σm-σn)^2)/R
end
