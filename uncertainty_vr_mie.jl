using Pkg

Pkg.activate("../..")
using Clapeyron, DataFrames, CSV, Distributions, Random
using PyCall
import PyPlot; const plt=PyPlot
corner = pyimport("corner")

function saturation_p_rhol(model::EoSModel,T)
    sat = saturation_pressure(model,T)
    return sat[1], 1/sat[2]
end

Nexp = 10

species = [
"decane"]

for s in species
    println(s)

    system = GERG2008([s])

    crit_exp = crit_pure(system)
    Tc = crit_exp[1]
    T_exp = range(0.3*Tc,0.9*Tc,length=Nexp)
    T_exp = [T_exp[i] for i in 1:Nexp]

    sat_exp = saturation_pressure.(system,T_exp)
    Hvap_exp = enthalpy_vap.(system,T_exp)

    psat_exp = [sat_exp[i][1] for i in 1:Nexp]
    rholsat_exp = [1/sat_exp[i][2] for i in 1:Nexp]

    a=Any["saturation_p_rhol","T"]
    append!(a,T_exp)
    b=Any["","out_pre"]
    append!(b,psat_exp)
    d=Any["","out_rho"]
    append!(d,rholsat_exp)

    df_sat = DataFrame(A=a, B=b, D=d)

    CSV.write("saturation_p_rho.csv",df_sat)

    model = SAFTVRMie([s])

    toestimate = [
        Dict(
            :param => :segment,
            :lower => model.params.segment.values[1]*0.9,
            :upper => model.params.segment.values[1]*1.1,
            :guess => model.params.segment.values[1]
        ),
        Dict(
            :param => :epsilon,
            :lower => model.params.epsilon.values[1]*0.9,
            :upper => model.params.epsilon.values[1]*1.1,
            :guess => model.params.epsilon.values[1]
        ),
        Dict(
            :param => :sigma,
            :factor => 1e-10,
            :lower => model.params.sigma.values[1]*0.9/1e-10,
            :upper => model.params.sigma.values[1]*1.1/1e-10,
            :guess => model.params.sigma.values[1]/1e-10
        ),
        Dict(
            :param => :lambda_r,
            :lower => model.params.lambda_r.values[1]*0.7,
            :upper => model.params.lambda_r.values[1]*1.5,
            :guess => model.params.lambda_r.values[1]
        ),
        Dict(
            :param => :lambda_a,
            :lower => model.params.lambda_a.values[1]*0.9,
            :upper => model.params.lambda_a.values[1]*1.2,
            :guess => model.params.lambda_a.values[1]
        )
    ]

    e = Estimation(model,toestimate,["saturation_p_rho.csv"])

    optimize!(e,Clapeyron.Metaheuristics.ECA())

    Nₛ = 1000
    Nᵢ = 1000

    θ₀ = [e.model.params.segment.values[1],
        e.model.params.epsilon.values[1],
        e.model.params.sigma.values[1]/1e-10,
        e.model.params.lambda_r.values[1],
        e.model.params.lambda_a.values[1]]
    σₚ = [0.08,0.08,0.08,0.1,0.1].*θ₀

    dₙ = Normal.(θ₀,σₚ)

    θᵢ = ones(Nᵢ,Nₛ,5)
    θᵢ[1,:,:] = ones(Nₛ,5).*θ₀'

    θ₂ = ones(Nₛ,5)

    dₛ = Chisq(Nexp-1)

    Δy₀ = [Clapeyron.obj_fun(e,θᵢ[1,i,:]) for i in 1:Nₛ]

    for i in 1:Nᵢ-1
        θᵢ[i+1,:,:] = deepcopy(θᵢ[i,:,:])

        θ₂[:,1] = rand(dₙ[1],Nₛ)
        θ₂[:,2] = rand(dₙ[2],Nₛ)
        θ₂[:,3] = rand(dₙ[3],Nₛ)
        θ₂[:,4] = rand(dₙ[4],Nₛ)
        θ₂[:,5] = rand(dₙ[5],Nₛ)

        Δy₁ = [Clapeyron.obj_fun(e,θ₂[i,:]) for i in 1:Nₛ]
        PDF = exp.(-Δy₁.+Δy₀)
        accept = ((PDF.>1) .| (PDF .>= rand(Nₛ)))
        θᵢ[i+1,accept,:] = deepcopy(θ₂[accept,:])
        Δy₀[accept] = Δy₁[accept]

        println("Iter: "*string(i)*", Fraction accepted: "*string(sum(accept)/length(accept)))
    end

    θᵢ = reshape(θᵢ,:,5)

    figure = corner.corner(
        θᵢ[100*1000:end,:],
        labels=[
            "m",
            "ϵ",
            "σ",
            "λᵣ",
            "λₐ"   
        ],
        quantiles=[0.025, 0.5, 0.975],
        show_titles=true,
        truths=[model.params.segment.values[1],
                model.params.epsilon.values[1],
                model.params.sigma.values[1]/1e-10,
                model.params.lambda_r.values[1],
                model.params.lambda_a.values[1]],
        truth_color="red"
    )

    figure.savefig("corner_vrmie_"*s*".pdf")

    figure.clf()
end