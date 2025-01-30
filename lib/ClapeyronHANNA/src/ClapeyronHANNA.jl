module ClapeyronHANNA

# using Reexport
using Clapeyron
const C = Clapeyron
using Clapeyron: EoSParam, ActivityModel, SingleParam, EoSVectorParam
using Clapeyron: getparams, init_puremodel, Rgas
import Clapeyron: default_locations
using Flux, Transformers.HuggingFace, Transformers.TextEncoders
using CSV, JSON3, JLD2
using LinearAlgebra

include("HANNA.jl")

const DB_PATH = normpath(Base.pkgdir(ClapeyronHANNA),"data")

# Constants
const N_emb_chembert = 384
const N_nodes = 96

"""
    load_chembert()

Load ChemBERTa model from HuggingFace (`DeepChem/ChemBERTa-77M-MTR`).
"""
function load_chembert(;name="DeepChem/ChemBERTa-77M-MTR", max_length=512, download=false)
    if download
        config = load_config(name)
        myconfig = HuggingFace.HGFConfig(config; max_length, hidden_act="gelu_noapprox")
        return load_model(name; config=myconfig), load_custom_tokenizer(max_length=max_length)
    else
        chembert = load("$DB_PATH/ChemBERTa/ChemBERTa-77M-MTR.jld2")
        return chembert["model"], load_custom_tokenizer(max_length=max_length)
    end
end

function load_custom_tokenizer(;max_length)
    vocab = JSON3.read(read("$DB_PATH/ChemBERTa/vocab.json",String))
    startsym = "[CLS]"
    endsym = "[SEP]"
    unksym = "[UNK]"
    padsym = "[PAD]"
    return TransformerTextEncoder(split_smiles, String.(keys(vocab)); startsym, endsym, unksym, padsym=padsym, trunc=max_length)
end

function split_smiles(smiles::String)
    return [match.match for match in eachmatch(r"\[(.*?)\]|Br|Cl|.", smiles)]
end

function load_scalers()
    data = CSV.File("$DB_PATH/scaler/scaler.csv";header=1,types=Float64,delim=',') |> CSV.Tables.matrix
    (u,s,k) = [data[:,i][:] for i in 1:3]
    T_scaler = TScaler(u[1],s[1]/k[1])
    emb_scaler = EmbeddedScaler(u,s,k)
    return T_scaler, emb_scaler
end

struct TScaler <: Function
    u::Float64
    sk::Float64
end

function (t::TScaler)(x)
    (x - t.u) / t.sk
end
 
struct EmbeddedScaler{D} <: Function
    u::D
    s::D
    k::D
end

function (t::EmbeddedScaler)(x)
    (x .- @view(t.u[2:end])) ./  @view(t.s[2:end]) .*  @view(t.k[2:end])
end


silu(x) = @. x/(1+exp(-x))

# Cosine similarity
function cosine_similarity(x1,x2;eps=1e-8)
    ∑x1 = sqrt(dot(x1,x1))
    ∑x2 = sqrt(dot(x2,x2))
    return dot(x1,x2)/(max(∑x1,eps*one(∑x1))*max(∑x2,eps*one(∑x2)))
end

end #module
