using Pkg
using Clapeyron, Metaheuristics
using DataFrames, CSV
using Distributions, Random
using PyCall
import PyPlot; const plt=PyPlot
corner = pyimport("corner")

function saturation_p(model::EoSModel,T)
    sat = saturation_pressure(model,T)
    return sat[1]
end

function saturation_rhol(model::EoSModel,T)
    sat = saturation_pressure(model,T)
    return 1/sat[2]
end

Nexp = 30

species = ["methane"]

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

    a=Any["[method=saturation_p]","T"]
    append!(a,T_exp)
    b=Any["","out_pre"]
    append!(b,psat_exp)

    df_psat = DataFrame(A=a, B=b)

    CSV.write("saturation_pressure.csv",df_psat)

    c=Any["[method=saturation_rhol]","T"]
    append!(c,T_exp)
    d=Any["","out_rho"]
    append!(d,rholsat_exp)

    df_rhosat = DataFrame(C=c, D=d)

    CSV.write("saturation_liquid_density.csv",df_rhosat)

    model = PCSAFT([s])

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
        )
    ]

    e,objective,initial,upper,lower = Estimation(model,toestimate,["saturation_liquid_density.csv","saturation_pressure.csv"])

    bounds = [lower upper]
    result = optimize(objective, bounds, ECA())
    params = minimizer(result)

    model = return_model(e,model,params)

    Nₛ = 1000
    Nᵢ = 1000

    θ₀ = [model.params.segment.values[1],model.params.epsilon.values[1],model.params.sigma.values[1]/1e-10]
    println(θ₀)
    σₚ = 0.08*θ₀

    dₙ = Normal.(θ₀,σₚ)

    θᵢ = ones(Nᵢ,Nₛ,3)
    θᵢ[1,:,:] = ones(Nₛ,3).*θ₀'

    θ₂ = ones(Nₛ,3)

    dₛ = Chisq(Nexp-1)

    Δy₀ = [objective(θᵢ[1,i,:]) for i in 1:Nₛ]

    for i in 1:Nᵢ-1
        θᵢ[i+1,:,:] = deepcopy(θᵢ[i,:,:])

        θ₂[:,1] = rand(dₙ[1],Nₛ)
        θ₂[:,2] = rand(dₙ[2],Nₛ)
        θ₂[:,3] = rand(dₙ[3],Nₛ)

        Δy₁ = [objective(θ₂[i,:]) for i in 1:Nₛ]
        PDF = exp.(-Δy₁.+Δy₀)
        accept = ((PDF.>1) .| (PDF .>= rand(Nₛ)))
        θᵢ[i+1,accept,:] = deepcopy(θ₂[accept,:])
        Δy₀[accept] = Δy₁[accept]
        println(sum(accept)/1000)
    end

    θᵢ = reshape(θᵢ,:,3)

    figure = corner.corner(
        θᵢ[100*1000:end,:],
        labels=[
            "m",
            "ϵ",
            "σ",
        ],
        quantiles=[0.025, 0.5, 0.975],
        show_titles=true,
        truths=[model.params.segment.values[1],model.params.epsilon.values[1],model.params.sigma.values[1]/1e-10],
        truth_color="red"
    )

    figure.savefig("corner_"*s*".pdf")

    figure.clf()
end