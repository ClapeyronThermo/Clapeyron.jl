using CSV, DataFrames
struct SAFTParam
    segment::Dict
    sigma::Dict
    epsilon::Dict
    lambdaA::Dict
    lambdaR::Dict
    k::Dict # keys are tuples of every pair
end

abstract type SAFT end

abstract type SAFTFamily <: SAFT end

struct system <: SAFTFamily; components; parameters::SAFTParam end

function System(species::AbstractArray,method::String = "PCSAFT")
    segment = Dict()
    sigma = Dict()
    epsilon = Dict()
    lambdaA = Dict()
    lambdaR = Dict()

    k = Dict()

    f = CSV.Rows("database/data_" * method * ".csv"; header=3)
    found_row = []
    i = 1
    for s in species
        count = 0
        for row in f
            count+=1
            if row.Compound == s
                push!(found_row, count)
                break
            end
        end
        dF = CSV.read("database/data_" * method * ".csv"; header=3, datarow=found_row[1]+3, limit=1)
        segment[i] = dF[1,:m]
        sigma[i,i] = dF[1,:sigma]*1e-10
        epsilon[i,i] = dF[1,:epsilon]
        if method == "SAFTVRMie"
            lambdaA[i,i] = dF[1,:lambdaA]
            lambdaR[i,i] = dF[1,:lambdaR]
        end
        k[i,i] = 0
        i += 1
    end

    model = system([1], SAFTParam(segment, sigma,epsilon,lambdaA,lambdaR, k)) # changed comp to [1, 2]

    include("Ideal.jl")
    include(method*".jl")
    include("Methods.jl")

    EoS(z,v,T)    = N_A*k_B*sum(z[i] for i in model.components)*T*(a_ideal(model,z,v,T)+a_res(model,z,v,T))
    return model,EoS
end
