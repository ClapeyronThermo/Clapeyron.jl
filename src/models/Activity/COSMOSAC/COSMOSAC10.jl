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

export COSMOSAC10

default_locations(::Type{COSMOSAC10}) = ["Activity/COSMOSAC/COSMOSAC10_like.csv"]

"""
    COSMOSAC10(components;
    puremodel = PR,
    userlocations = String[],
    pure_userlocations = String[],
    verbose = false)

## Input parameters:
- `Pnhb` :Single Parameter{String}
- `POH` :Single Parameter{String}
- `POT` :Single Parameter{String}
- `V`: Single Parameter{Float64}
- `A`: Single Parameter{Float64}

## Description
An activity coefficient model using molecular solvation based on the COSMO-RS method. Sigma profiles are now split by non-hydrogen bonding, hydrogen acceptor and hydrogen donor.

## References
1. Klamt, A. (1995). Conductor-like screening model for real solvents: A new approach to the quantitative calculation of solvation phenomena. Journal of Physical Chemistry, 99(7), 2224–2235. [doi:10.1021/j100007a062](https://doi.org/10.1021/j100007a062)
2. Lin, S-T. & Sandler, S.I. (2002). A priori phase equilibrium prediction from a segment contribution solvation model. Industrial & Engineering Chemistry Research, 41(5), 899–913. [doi:10.1021/ie001047w](https://doi.org/10.1021/ie001047w)
3. Hsieh, C-H., Sandler, S.I., & Lin, S-T. (2010). Improvements of COSMO-SAC for vapor–liquid and liquid–liquid equilibrium predictions. Fluid Phase Equilibria, 297(1), 90-97. [doi:10.1016/j.fluid.2010.06.011](https://doi.org/10.1016/j.fluid.2010.06.011)
"""
COSMOSAC10

function COSMOSAC10(components;
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
        Pnhb = [zeros(51) for i in 1:length(components)]
        POH = [zeros(51) for i in 1:length(components)]
        POT = [zeros(51) for i in 1:length(components)]
        for i in 1:length(components)
            id = cas(formatted_components[i])
            ids = CAS.==uppercase(id[1])
            dbname = INCHIKEY[ids]
            file = String(take!(Downloads.download("https://raw.githubusercontent.com/usnistgov/COSMOSAC/master/profiles/UD/sigma3/"*dbname[1]*".sigma", IOBuffer())))
            lines = split(file,r"\n")
            meta = lines[1][9:end]
            json = JSON3.read(meta)
            A[i] = json["area [A^2]"]
            V[i] = json["volume [A^3]"]
            Pnhb[i] = [parse(Float64,split(lines[i]," ")[2]) for i in 4:54]
            POH[i] = [parse(Float64,split(lines[i]," ")[2]) for i in 55:105]
            POT[i] = [parse(Float64,split(lines[i]," ")[2]) for i in 106:156]
        end
        A = SingleParam("A",formatted_components,A)
        V = SingleParam("V",formatted_components,V)
        Pnhb = SingleParam("Pnhb",formatted_components,Pnhb)
        POH = SingleParam("POH",formatted_components,POH)
        POT = SingleParam("POT",formatted_components,POT)
    else
        params = getparams(formatted_components, default_locations(COSMOSAC10); userlocations = userlocations, ignore_missing_singleparams=["Pnhb","POH","POT","A","V"], verbose = verbose)
        Pnhb  = COSMO_parse_Pi(params["Pnhb"])
        POH  = COSMO_parse_Pi(params["POH"])
        POT  = COSMO_parse_Pi(params["POT"])
        A  = params["A"]
        V  = params["V"]
    end

    _puremodel = init_puremodel(puremodel,components,pure_userlocations,verbose)
    packagedparams = COSMOSAC10Param(Pnhb,POH,POT,V,A)
    references = ["10.1021/acs.jctc.9b01016","10.1021/acs.iecr.7b01360","10.1021/j100007a062"]
    model = COSMOSAC10(formatted_components,packagedparams,_puremodel,1e-12,references)
    return model
end

function Γ_as_view(Γ,l1 = length(Γ) ÷ 3)
    Γnhb = viewn(Γ, 51, 1)
    ΓOH = viewn(Γ, 51, 2)
    ΓOT = viewn(Γ, 51, 3)
    return Γnhb, ΓOH, ΓOT
end

function excess_g_res(model::COSMOSAC10Model,V,T,z)
    lnγ = @f(lnγ_res)
    sum(z[i]*R̄*T*lnγ[i] for i ∈ @comps)
end

function lnγ_res(model::COSMOSAC10Model,V,T,z)   
    A = model.params.A.values
    Ā = dot(z,A)
    PS = zeros(typeof(Ā),51*3)
    PSᵢ = PS
    PSnhb, PSOH, PSOT = Γ_as_view(PS,51)
    Pnhb  = model.params.Pnhb.values
    POH = model.params.POH.values
    POT = model.params.POT.values
    @inbounds @simd for v in 1:51
        PS_nhbᵢᵥ = zero(eltype(PS))
        PS_OHᵢᵥ = zero(eltype(PS))
        PS_OTᵢᵥ = zero(eltype(PS))
        for i in @comps
            zᵢ = z[i]
            PS_nhbᵢᵥ += Pnhb[i][v]*zᵢ
            PS_OHᵢᵥ += POH[i][v]*zᵢ
            PS_OTᵢᵥ += POT[i][v]*zᵢ
        end
        PSnhb[v] = PS_nhbᵢᵥ
        PSOH[v] = PS_OHᵢᵥ
        PSOT[v] = PS_OTᵢᵥ
    end
    PS ./= Ā
    #n = A ./ aeff
    lnΓS = @f(lnΓ,PS)
    (lnΓSnhb, lnΓSOH, lnΓSOT)= lnΓS
    lnΓi = Vector{typeof(lnΓS)}(undef,length(model))
    
    PSnhbᵢ, PSOHᵢ, PSOTᵢ = Γ_as_view(PSᵢ,51)
    #lnΓi = [@f(lnΓ,Pnhb[i]./A[i],POH[i]./A[i],POT[i]./A[i]) for i ∈ @comps]
    for i in @comps
        Aᵢ = A[i]
        PSnhbᵢ .= Pnhb[i] ./ Aᵢ
        PSOHᵢ .= POH[i] ./ Aᵢ
        PSOTᵢ .= POT[i] ./ Aᵢ
        lnΓi[i] = @f(lnΓ,PSᵢ)
    end
    aeff = 7.5
    aeff⁻¹ = 1/7.5
    lnγ_res = zeros(eltype(lnΓSnhb),length(model))
    for i in @comps
        #nᵢ = A[i]/aeff
        #Aᵢ⁻¹ = 1/A[i]
        lnγ_resᵢ = zero(eltype(lnγ_res))
        lnΓSnhbᵢ, lnΓSOHᵢ, lnΓSOTᵢ = lnΓi[i]
        Pnhbᵢ, POHᵢ, POTᵢ = Pnhb[i], POH[i], POT[i]
        for v in 1:51
            lnγ_resᵢ += aeff⁻¹*Pnhbᵢ[v]*(lnΓSnhb[v] - lnΓSnhbᵢ[v])
            lnγ_resᵢ += aeff⁻¹*POHᵢ[v]*(lnΓSOH[v] - lnΓSOHᵢ[v])
            lnγ_resᵢ += aeff⁻¹*POTᵢ[v]*(lnΓSOT[v] - lnΓSOTᵢ[v])
        end
        lnγ_res[i] = lnγ_resᵢ
    end
    return lnγ_res
    #=
    
    lnγ_res = [A[i]/aeff*(sum(Pnhb[i][v]/A[i]*(lnΓSnhb[v]-lnΓi[i][1][v]) for v ∈ 1:51)
                      +sum(POH[i][v]/A[i]*(lnΓSOH[v]-lnΓi[i][2][v]) for v ∈ 1:51)
                      +sum(POT[i][v]/A[i]*(lnΓSOT[v]-lnΓi[i][3][v]) for v ∈ 1:51)) for i ∈ @comps]
    
    return lnγ_res =#
end

function lnΓ(model::COSMOSAC10Model,V,T,z,P)
    Γ = _lnΓ(model,V,T,z,P)
    return Γ_as_view(Γ,51)
end

function PΔW(model::COSMOSAC10Model,V,T,z,P,nonzeros = 1:length(P))
    σ  = -0.025:0.001:0.025
    Tinv = 1/T
    TYPE = @f(Base.promote_eltype,P,Tinv)
    PW = zeros(TYPE,length(nonzeros),length(nonzeros))
    idx(i) = div(i-1,51) + 1
    #σ(i) = _σ[rem(i + 50,51) + 1]
    idxx(i) = rem(i + 50,51) + 1
    i_nonzero = 0
    @inbounds for ii in 1:153
        i = idxx(ii)
        iszero(P[ii]) && continue
        i_nonzero +=1
        v_nonzero = 0
        for vv in 1:153
            iszero(P[vv]) && continue
            v_nonzero +=1
            v = idxx(vv)
            PW[i_nonzero,v_nonzero] = P[vv]*exp(-ΔW(σ[i],σ[v],idx(vv),idx(ii),T)*Tinv)
        end
    end
    return PW
end

function ΔW(σm,σn,t,s,T)
    ces  = 6525.69+1.4859e8/T^2
    chb = COSMOSAC10_ΔW_data[t,s] * (σm*σn<0)
    R  = 0.001987
    return (ces*(σm+σn)^2-chb*(σm-σn)^2)/R
end

const COSMOSAC10_ΔW_data =@SMatrix [0.0 0.0 0.0;0.0 4013.78 0.0;0.0 3016.43 932.31]

#fcosmo(system::COSMOSAC10) = Clapeyron.activity_coefficient(system,1e5, 333.15,[0.5,0.5])[1] - 1.4015660588643404
