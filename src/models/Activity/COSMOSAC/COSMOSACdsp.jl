struct COSMOSACdspParam <: EoSParam
    Pnhb::SingleParam{Vector{Float64}}
    POH::SingleParam{Vector{Float64}}
    POT::SingleParam{Vector{Float64}}
    epsilon::SingleParam{Float64}
    V::SingleParam{Float64}
    A::SingleParam{Float64}
    water::SingleParam{Int64}
    COOH::SingleParam{Int64}
    hb_acc::SingleParam{Int64}
    hb_don::SingleParam{Int64}
end

abstract type COSMOSACdspModel <: COSMOSAC10Model end

struct COSMOSACdsp{c<:EoSModel} <: COSMOSACdspModel
    components::Array{String,1}
    params::COSMOSACdspParam
    puremodel::EoSVectorParam{c}
    absolutetolerance::Float64
    references::Array{String,1}
end

export COSMOSACdsp

default_locations(::Type{COSMOSACdsp}) = ["Activity/COSMOSAC/COSMOSAC10_like.csv","Activity/COSMOSAC/COSMOSACdsp_like.csv"]

function COSMOSACdsp(components;
    puremodel = PR,
    userlocations = String[],
    pure_userlocations = String[],
    verbose=false)

    formatted_components = format_components(components)
    params = getparams(formatted_components, default_locations(COSMOSACdsp); userlocations=userlocations, ignore_missing_singleparams=["Pnhb","POH","POT","A","V","epsilon","water","COOH","hb_acc","hb_don"], verbose=verbose)
    Pnhb  = COSMO_parse_Pi(params["Pnhb"])
    POH  = COSMO_parse_Pi(params["POH"])
    POT  = COSMO_parse_Pi(params["POT"])
    A  = params["A"]
    V  = params["V"]
    epsilon  = params["epsilon"]
    water = params["water"]
    COOH = params["COOH"]
    hb_acc = params["hb_acc"]
    hb_don = params["hb_don"]

    if any(V.values .==0)
        CAS, INCHIKEY = get_cosmo_comps()
        for i in 1:length(components)
            if V.values[i]==0
                id = cas(components[i])
                ids = CAS.==uppercase(id[1])
                dbname = INCHIKEY[ids]
                file = String(take!(Downloads.download("https://raw.githubusercontent.com/usnistgov/COSMOSAC/master/profiles/UD/sigma/"*dbname[1]*".sigma", IOBuffer())))
                lines = split(file,r"\n")
                meta = lines[1][9:end]
                json = JSON3.read(meta)
                A.values[i] = json["area [A^2]"]
                A.ismissingvalues[i] = false
                V.values[i] = json["volume [A^3]"]
                V.ismissingvalues[i] = false
                Pnhb.values[i] = [parse(Float64,split(lines[i]," ")[2]) for i in 4:54]
                Pnhb.ismissingvalues[i] = false
                POH.values[i] = [parse(Float64,split(lines[i]," ")[2]) for i in 55:105]
                POH.ismissingvalues[i] = false
                POT.values[i] = [parse(Float64,split(lines[i]," ")[2]) for i in 106:156]
                POT.ismissingvalues[i] = false
                epsilon.values[i] = json["disp. e/kB [K]"]
                epsilon.ismissingvalues[i] = false
                if json["disp. flag"] == "NHB"
                    COOH.values[i] = 0
                    COOH.ismissingvalues[i] = false
                    hb_acc.values[i] = 0
                    hb_acc.ismissingvalues[i] = false
                    hb_don.values[i] = 0
                    hb_don.ismissingvalues[i] = false
                    water.values[i] = 0
                    water.ismissingvalues[i] = false
                elseif json["disp. flag"] == "COOH"
                    COOH.values[i] = 1
                    COOH.ismissingvalues[i] = false
                    hb_acc.values[i] = 0
                    hb_acc.ismissingvalues[i] = false
                    hb_don.values[i] = 0
                    hb_don.ismissingvalues[i] = false
                    water.values[i] = 0
                    water.ismissingvalues[i] = false
                elseif json["disp. flag"] == "HB-DONOR-ACCEPTOR"
                    COOH.values[i] = 0
                    COOH.ismissingvalues[i] = false
                    hb_acc.values[i] = 1
                    hb_acc.ismissingvalues[i] = false
                    hb_don.values[i] = 1
                    hb_don.ismissingvalues[i] = false
                    water.values[i] = 0
                    water.ismissingvalues[i] = false
                elseif json["disp. flag"] == "HB-ACCEPTOR"
                    COOH.values[i] = 0
                    COOH.ismissingvalues[i] = false
                    hb_acc.values[i] = 1
                    hb_acc.ismissingvalues[i] = false
                    hb_don.values[i] = 0
                    hb_don.ismissingvalues[i] = false
                    water.values[i] = 0
                    water.ismissingvalues[i] = false
                elseif json["disp. flag"] == "HB-DONOR"
                    COOH.values[i] = 0
                    COOH.ismissingvalues[i] = false
                    hb_acc.values[i] = 0
                    hb_acc.ismissingvalues[i] = false
                    hb_don.values[i] = 1
                    hb_don.ismissingvalues[i] = false
                    water.values[i] = 0
                    water.ismissingvalues[i] = false
                end
            end
        end
    end

    _puremodel = init_puremodel(puremodel,components,pure_userlocations,verbose)
    packagedparams = COSMOSACdspParam(Pnhb,POH,POT,epsilon,V,A,water,COOH,hb_acc,hb_don)
    references = String[]
    model = COSMOSACdsp(formatted_components,packagedparams,_puremodel,1e-12,references)
    return model
end

function activity_coefficient(model::COSMOSACdspModel,V,T,z)
    return exp.(@f(lnγ_comb) .+@f(lnγ_res).+@f(lnγ_dsp))
end

function excess_g_res(model::COSMOSACdspModel,V,T,z)
    lnγ = @f(lnγ_res)
    sum(z[i]*R̄*T*lnγ[i] for i ∈ @comps)
end

function lnγ_dsp(model::COSMOSACdspModel,V,T,z)
    x = z./sum(z)
    if model.params.water.values[1]==1 && model.params.hb_acc.values[2]==1 && model.params.hb_don.values[2]==0
        w = -0.27027
    elseif model.params.water.values[2]==1 && model.params.hb_acc.values[1]==1 && model.params.hb_don.values[1]==0
        w = -0.27027
    elseif model.params.COOH.values[1]==1 && (model.params.hb_acc.values[2]==1 | model.params.hb_don.values[2]==1)
        w = -0.27027
    elseif model.params.COOH.values[1]==1 && model.params.water.values[2]==1
        w = -0.27027
    elseif model.params.COOH.values[2]==1 && model.params.water.values[1]==1
        w = -0.27027
    else
        w = 0.27027
    end

    ϵ = model.params.epsilon.values

    A = w*(0.5*(ϵ[1]+ϵ[2])-√(ϵ[1]*ϵ[2]))
    return A*(1 .-x).^2
end
