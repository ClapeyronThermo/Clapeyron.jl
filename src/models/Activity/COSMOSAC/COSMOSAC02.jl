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
    verbose = false,
    reference_state = nothing)

## Input parameters:
- `Pi` :Single Parameter{String} 
- `V`: Single Parameter{Float64}
- `A`: Single Parameter{Float64}

## Description
An activity coefficient model using molecular solvation based on the COSMO-RS method.

## References
1. Klamt, A. (1995). Conductor-like screening model for real solvents: A new approach to the quantitative calculation of solvation phenomena. Journal of Physical Chemistry, 99(7), 2224–2235. [doi:10.1021/j100007a062](https://doi.org/10.1021/j100007a062)
2. Lin, S-T. & Sandler, S.I. (2002). A priori phase equilibrium prediction from a segment contribution solvation model. Industrial & Engineering Chemistry Research, 41(5), 899–913. [doi:10.1021/ie001047w](https://doi.org/10.1021/ie001047w)
"""
COSMOSAC02

default_locations(::Type{COSMOSAC02}) = ["Activity/COSMOSAC/COSMOSAC02_like.csv"]

function COSMOSAC02(components;
    puremodel = PR,
    userlocations = String[],
    pure_userlocations = String[],
    use_nist_database = false,
    verbose = false,
    reference_state = nothing)
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
    references = String["10.1021/ie001047w","10.1021/acs.jctc.9b01016","10.1021/acs.iecr.7b01360","10.1021/j100007a062"]
    model = COSMOSAC02(formatted_components,packagedparams,_puremodel,1e-12,references)
    set_reference_state!(model,reference_state,verbose = verbose)
    return model
end

function activity_coefficient(model::COSMOSAC02Model,V,T,z)
    return exp.(@f(lnγ_comb) .+ @f(lnγ_res))
end

function lnγ_impl!(lnγ,model::COSMOSAC02Model,V,T,z)
    lnγ .= 0
    lnγ .+= @f(lnγ_res)
    lnγ .+= @f(lnγ_comb)
    return lnγ
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

lnΓ(model::COSMOSAC02Model,V,T,z,P) = @f(_lnΓ,P)

function _lnΓ(model::COSMOSAC02Model,V,T,z,P)
    #=
    we perform three strategies to solve the system
    1. we preallocate P(sigma)*W(sigma)
    2. we remove P(sigma) = 0 entries and iterate on the reduce matrix.
    3. we use the last calculated  Γ[i] in the fixpoint iteration.
    =#
    _TYPE = @f(Base.promote_eltype)
    #strategy 2
    nonzeros = findall(!iszero,P)
    l1 = length(nonzeros)
    Γ0 = ones(_TYPE,l1)
    #strategy 1
    PW = @f(PΔW,P,nonzeros)
    
    function f!(Γnew,Γold)
        i_solved = 0
        for ii in 1:length(Γold)
            _res = zero(eltype(Γnew))
            #strategy 3
            @inbounds for vv in 1:i_solved
                _res += Γnew[vv]*PW[ii,vv]
            end
            @inbounds for vv in (i_solved+1):length(Γold)
                _res += Γold[vv]*PW[ii,vv]
            end
            Γnew[ii] = 1/_res
            i_solved += 1
        end
        #this would be the code without strategy 3
        #mul!(Γnew,PW,Γold)
        #Γnew .= 1 ./ Γnew
        return Γnew
    end
    Γres = Solvers.fixpoint(f!,Γ0,Solvers.SSFixPoint(dampingfactor = 0.5,lognorm = true,normorder = 1),max_iters = 500*length(model),atol = model.absolutetolerance,rtol = 0.0)
    Γres .= log.(Γres)
    #restore original array size.
    resize!(Γ0,length(P))
    Γ0 .= 0.0
    for (k,i) in pairs(nonzeros)
        Γ0[i] = Γres[k]
    end
    return Γ0
end

function PΔW(model::COSMOSAC02Model,V,T,z,P,nonzeros = 1:length(P))
    _σ = -0.025:0.001:0.025
    σ  = @view _σ[nonzeros]
    Tinv = 1/T
    TYPE = @f(Base.promote_eltype,P,Tinv)
    l = length(nonzeros)
    _P = @view P[nonzeros]    
    PW = zeros(TYPE,(l,l))
    @inbounds for i in 1:l
        for j in 1:l
            PW[i,j] = _P[j]*exp(-ΔW(σ[j],σ[i])*Tinv) 
        end
    end
    return PW
end

function ΔW(σm,σn)
    σdon,σacc = minmax(σm,σn)
    σhb  = 0.0084
    chb  = 85580
    α    = 16466.72
    R    = 0.001987
    return (α/2*(σm+σn)^2+chb*max(0,σacc-σhb)*min(0,σdon+σhb))/R
end
