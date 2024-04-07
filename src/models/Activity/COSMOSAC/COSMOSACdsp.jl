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

"""
    COSMOSACdsp(components;
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
- `epsilon`: Single Parameter{Float64}
- `COOH`: Single Parameter{Float64}
- `water`: Single Parameter{Float64}
- `hb_acc`: Single Parameter{Float64}
- `hb_don`: Single Parameter{Float64}

## Description
An activity coefficient model using molecular solvation based on the COSMO-RS method. Sigma profiles are now split by non-hydrogen bonding, hydrogen acceptor and hydrogen donor. A dispersive correction is included.

## References
1. Hsieh, C-H., Lin, S-T. & Vrabec, J. (2014). Considering the dispersive interactions in the COSMO-SAC model for more accurate predictions of fluid phase behavior. Fluid Phase Equilibria, 367, 109-116. [doi:10.1016/j.fluid.2014.01.032](https://doi.org/10.1016/j.fluid.2014.01.032)
"""
COSMOSACdsp

function COSMOSACdsp(components;
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
            epsilon[i] = json["disp. e/kB [K]"]
            if json["disp. flag"] == "NHB"
                COOH[i] = 0
                hb_acc[i] = 0
                hb_don[i] = 0
                water[i] = 0
            elseif json["disp. flag"] == "COOH"
                COOH[i] = 1
                hb_acc[i] = 0
                hb_don[i] = 0
                water[i] = 0
            elseif json["disp. flag"] == "HB-DONOR-ACCEPTOR"
                COOH[i] = 0
                hb_acc[i] = 1
                hb_don[i] = 1
                water[i] = 0
            elseif json["disp. flag"] == "HB-ACCEPTOR"
                COOH[i] = 0
                hb_acc[i] = 1
                hb_don[i] = 0
                water[i] = 0
            elseif json["disp. flag"] == "HB-DONOR"
                COOH[i] = 0
                hb_acc[i] = 0
                hb_don[i] = 1
                water[i] = 0
            elseif json["disp. flag"] == "WATER"
                COOH[i] = 0
                hb_acc[i] = 0
                hb_don[i] = 0
                water[i] = 1
            end
        end
        A = SingleParam("A",formatted_components,A)
        V = SingleParam("V",formatted_components,V)
        Pnhb = SingleParam("Pnhb",formatted_components,Pnhb)
        POH = SingleParam("POH",formatted_components,POH)
        POT = SingleParam("POT",formatted_components,POT)
        epsilon = SingleParam("epsilon",formatted_components,epsilon)
        COOH = SingleParam("COOH",formatted_components,COOH)
        hb_acc = SingleParam("hb_acc",formatted_components,hb_acc)
        hb_don = SingleParam("hb_don",formatted_components,hb_don)
        water = SingleParam("water",formatted_components,water)
    else
        params = getparams(formatted_components, default_locations(COSMOSACdsp); userlocations = userlocations, ignore_missing_singleparams=["Pnhb","POH","POT","A","V","epsilon","water","COOH","hb_acc","hb_don"], verbose = verbose)
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
    end

    _puremodel = init_puremodel(puremodel,components,pure_userlocations,verbose)
    packagedparams = COSMOSACdspParam(Pnhb,POH,POT,epsilon,V,A,water,COOH,hb_acc,hb_don)
    references = ["10.1021/acs.jctc.9b01016","10.1021/acs.iecr.7b01360"]
    model = COSMOSACdsp(formatted_components,packagedparams,_puremodel,1e-12,references)
    return model
end

function activity_coefficient(model::COSMOSACdspModel,V,T,z)
    return exp.(@f(lnγ_comb) .+@f(lnγ_res).+@f(lnγ_dsp))
end

function excess_g_res(model::COSMOSACdspModel,V,T,z)
    lnγ = @f(lnγ_res) .+ @f(lnγ_dsp)
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
