using Clapeyron, DataFrames, CSV
using PyCall
import PyPlot; const plt = PyPlot

Nexp = 20

GERG = [
        "water"
        ]


for species in GERG
    system = GERG2008([species])

    crit_exp = crit_pure(system)
    Tc = crit_exp[1]
    T_exp = range(0.7*Tc,0.9*Tc,length=Nexp)
    T_exp = [T_exp[i] for i in 1:Nexp]

    sat_exp = saturation_pressure.(system,T_exp)
    Hvap_exp = enthalpy_vap.(system,T_exp)

    psat_exp = [sat_exp[i][1] for i in 1:Nexp]
    rholsat_exp = [1/sat_exp[i][2] for i in 1:Nexp]

    a=Any["saturation_p","T"]
    append!(a,T_exp)
    b=Any["","out_pre"]
    append!(b,psat_exp)

    df_psat = DataFrame(A=a, B=b)

    CSV.write("saturation_pressure.csv",df_psat)

    c=Any["saturation_rhol","T"]
    append!(c,T_exp)
    d=Any["","out_rho"]
    append!(d,rholsat_exp)

    df_rhosat = DataFrame(C=c, D=d)

    CSV.write("saturation_liquid_density.csv",df_rhosat)

    c=Any["enthalpy_vap","T"]
    append!(c,T_exp)
    d=Any["","out_Hvap"]
    append!(d,Hvap_exp)

    df_hvap = DataFrame(C=c, D=d)

    CSV.write("enthalpy_vap.csv",df_hvap)

    c=Any["crit_temp","out_crit"]
    append!(c,Tc)

    df_crit = DataFrame(C=c)

    CSV.write("crit_temp.csv",df_crit)

    model = PCSAFT([species])

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
            :param => :epsilon_assoc,
            :lower => model.params.epsilon_assoc.values.values[1]*0.9,
            :upper => model.params.epsilon_assoc.values.values[1]*1.1,
            :guess => model.params.epsilon_assoc.values.values[1]
        ),
        Dict(
            :param => :bondvol,
            :lower => model.params.bondvol.values.values[1]*0.9,
            :upper => model.params.bondvol.values.values[1]*1.1,
            :guess => model.params.bondvol.values.values[1]
        )
    ]

    e = Estimation(model,toestimate,["enthalpy_vap.csv"])

    # optimize!(e,Clapeyron.Metaheuristics.SA(x_initial=[e.toestimate.guess[1][1],e.toestimate.guess[2][1],e.toestimate.guess[3][1]],N=1000,tol_fun=1e-3))
    optimize!(e,Clapeyron.Metaheuristics.ECA())

    N = 1000

    mˢ = e.model.params.segment.values[1]
    ϵˢ = e.model.params.epsilon.values[1]
    σˢ = e.model.params.sigma.values[1]/1e-10
    Eˢ = e.model.params.epsilon_assoc.values.values[1]
    bˢ = e.model.params.bondvol.values.values[1]


    m = range(0.9*mˢ,1.1*mˢ,length=N)
    ϵ = range(0.9*ϵˢ,1.1*ϵˢ,length=N)
    σ = range(0.9*σˢ,1.1*σˢ,length=N)
    E = range(0.9*Eˢ,1.1*Eˢ,length=N)
    b = range(0.9*bˢ,1.1*bˢ,length=N)

    objm = zeros(N)
    objϵ = zeros(N)
    objσ = zeros(N)
    objE = zeros(N)
    objb = zeros(N)
    for i in 1:N
        objm[i] = Clapeyron.obj_fun(e,[m[i],ϵˢ,σˢ,Eˢ,bˢ])
        objϵ[i] = Clapeyron.obj_fun(e,[mˢ,ϵ[i],σˢ,Eˢ,bˢ])
        objσ[i] = Clapeyron.obj_fun(e,[mˢ,ϵˢ,σ[i],Eˢ,bˢ])
        objE[i] = Clapeyron.obj_fun(e,[mˢ,ϵˢ,σˢ,E[i],bˢ])
        objb[i] = Clapeyron.obj_fun(e,[mˢ,ϵˢ,σˢ,Eˢ,b[i]])

    end

    fig, axs = plt.subplots(5, 5)
    fig.subplots_adjust(wspace=0,hspace=0)
    fig.suptitle(species)
    axs[1,1].set_title("m")
    axs[1,1].scatter(m./mˢ.-1,objm,1,c=log10.(objm),vmax=1.5)
    axs[1,1].set_xlim(-0.1,0.1)
    axs[1,1].set_ylim(0,maximum(objm))
    axs[1,1].set_yticks([])
    axs[1,1].set_xticklabels([])
    axs[1,1].tick_params(axis="x",direction="in")
    axs[1,1].tick_params(axis="y",direction="in")

    axs[2,2].set_title("ϵ")
    axs[2,2].scatter(ϵ./ϵˢ.-1,objϵ,1,c=log10.(objϵ),vmax=1.5)
    axs[2,2].set_xlim(-0.1,0.1)
    axs[2,2].set_xticklabels([])
    axs[2,2].set_yticks([])
    axs[2,2].set_ylim(0,maximum(objϵ))
    axs[2,2].tick_params(axis="x",direction="in")
    axs[2,2].tick_params(axis="y",direction="in")

    axs[3,3].set_title("σ")
    axs[3,3].scatter(σ./σˢ.-1,objσ,1,c=log10.(objσ),vmax=1.5)
    axs[3,3].set_xlim(-0.1,0.1)
    axs[3,3].set_xticklabels([])
    axs[3,3].set_yticks([])
    axs[3,3].set_ylim(0,maximum(objσ))
    axs[3,3].tick_params(axis="x",direction="in")
    axs[3,3].tick_params(axis="y",direction="in")

    axs[4,4].set_title("ϵassoc")
    axs[4,4].scatter(σ./σˢ.-1,objσ,1,c=log10.(objσ),vmax=1.5)
    axs[4,4].set_xlim(-0.1,0.1)
    axs[4,4].set_xticklabels([])
    axs[4,4].set_yticks([])
    axs[4,4].set_ylim(0,maximum(objσ))
    axs[4,4].tick_params(axis="x",direction="in")
    axs[4,4].tick_params(axis="y",direction="in")

    axs[5,5].set_title("κ")
    axs[5,5].scatter(σ./σˢ.-1,objσ,1,c=log10.(objσ),vmax=1.5)
    axs[5,5].set_xlim(-0.1,0.1)
    axs[5,5].set_yticks([])
    axs[5,5].set_ylim(0,maximum(objσ))
    axs[5,5].set_xlabel("κ")
    axs[5,5].tick_params(axis="x",direction="in")
    axs[5,5].tick_params(axis="y",direction="in")
    axs[5,5].set_xticklabels([" ","0.0","0.1"])


    N = 100

    m = range(0.9*mˢ,1.1*mˢ,length=N).*ones(N,N)
    ϵ = range(0.9*ϵˢ,1.1*ϵˢ,length=N).*ones(N,N)
    σ = range(0.9*σˢ,1.1*σˢ,length=N).*ones(N,N)
    E = range(0.9*Eˢ,1.1*Eˢ,length=N).*ones(N,N)
    b = range(0.9*bˢ,1.1*bˢ,length=N).*ones(N,N)

    objmϵ = zeros(N,N)
    objmσ = zeros(N,N)
    objmE = zeros(N,N)
    objmb = zeros(N,N)
    objϵσ = zeros(N,N)
    objϵE = zeros(N,N)
    objϵb = zeros(N,N)
    objσE = zeros(N,N)
    objσb = zeros(N,N)
    objEb = zeros(N,N)

    for i in 1:N
        for j in 1:N
            objmϵ[i,j] = Clapeyron.obj_fun(e,[m[i,j],ϵ[j,i],σˢ,Eˢ,bˢ])
            objmσ[i,j] = Clapeyron.obj_fun(e,[m[i,j],ϵˢ,σ[j,i],Eˢ,bˢ])
            objmE[i,j] = Clapeyron.obj_fun(e,[m[i,j],ϵˢ,σˢ,E[j,i],bˢ])
            objmb[i,j] = Clapeyron.obj_fun(e,[m[i,j],ϵˢ,σˢ,Eˢ,b[j,i]])

            objϵσ[i,j] = Clapeyron.obj_fun(e,[mˢ,ϵ[i,j],σ[j,i],Eˢ,bˢ])
            objϵE[i,j] = Clapeyron.obj_fun(e,[mˢ,ϵ[i,j],σˢ,E[j,i],bˢ])
            objϵb[i,j] = Clapeyron.obj_fun(e,[mˢ,ϵ[i,j],σˢ,Eˢ,b[j,i]])

            objσE[i,j] = Clapeyron.obj_fun(e,[mˢ,ϵˢ,σ[i,j],E[j,i],bˢ])
            objσb[i,j] = Clapeyron.obj_fun(e,[mˢ,ϵˢ,σ[i,j],Eˢ,b[j,i]])

            objEb[i,j] = Clapeyron.obj_fun(e,[mˢ,ϵˢ,σˢ,E[i,j],b[j,i]])
        end
        println(i)
    end

    axs[2,1].pcolor(m./mˢ.-1,(ϵ./ϵˢ.-1)',log10.(objmϵ),vmax=1.5)
    axs[2,1].set_xticklabels([])
    axs[2,1].set_yticklabels([" "," ","0.0","0.1"])
    axs[2,1].set_ylabel("ϵ")
    axs[2,1].tick_params(axis="x",direction="in")
    axs[2,1].tick_params(axis="y",direction="in")

    axs[3,1].pcolor(m./mˢ.-1,(σ./σˢ.-1)',log10.(objmσ),vmax=1.5)
    axs[3,1].set_xticklabels([])
    axs[3,1].set_yticklabels([" "," ","0.0","0.1"])
    axs[3,1].set_ylabel("σ")
    axs[3,1].tick_params(axis="x",direction="in")
    axs[3,1].tick_params(axis="y",direction="in")

    axs[4,1].pcolor(m./mˢ.-1,(E./Eˢ.-1)',log10.(objmE),vmax=1.5)
    axs[4,1].set_xticklabels([])
    axs[4,1].set_yticklabels([" "," ","0.0","0.1"])
    axs[4,1].set_ylabel("ϵassoc")
    axs[4,1].tick_params(axis="x",direction="in")
    axs[4,1].tick_params(axis="y",direction="in")

    axs[5,1].pcolor(m./mˢ.-1,(b./bˢ.-1)',log10.(objmE),vmax=1.5)
    axs[5,1].set_xlabel("m")
    axs[5,1].set_ylabel("κ")
    axs[5,1].tick_params(axis="x",direction="in")
    axs[5,1].tick_params(axis="y",direction="in")

    axs[3,2].pcolor(ϵ./ϵˢ.-1,(σ./σˢ.-1)',log10.(objϵσ),vmax=1.5)
    axs[3,2].set_yticklabels([])
    axs[3,2].set_xticklabels([])
    axs[3,2].tick_params(axis="x",direction="in")
    axs[3,2].tick_params(axis="y",direction="in")

    axs[4,2].pcolor(ϵ./ϵˢ.-1,(E./Eˢ.-1)',log10.(objϵE),vmax=1.5)
    axs[4,2].set_yticklabels([])
    axs[4,2].set_xticklabels([])
    axs[4,2].tick_params(axis="x",direction="in")
    axs[4,2].tick_params(axis="y",direction="in")

    axs[5,2].pcolor(ϵ./ϵˢ.-1,(b./bˢ.-1)',log10.(objϵb),vmax=1.5)
    axs[5,2].set_yticklabels([])
    axs[5,2].set_xlabel("ϵ")
    axs[5,2].set_xticklabels([" "," ","0.0","0.1"])
    axs[5,2].tick_params(axis="x",direction="in")
    axs[5,2].tick_params(axis="y",direction="in")

    axs[4,3].pcolor(σ./σˢ.-1,(E./Eˢ.-1)',log10.(objσE),vmax=1.5)
    axs[4,3].set_yticklabels([])
    axs[4,3].set_xticklabels([])
    axs[4,3].tick_params(axis="x",direction="in")
    axs[4,3].tick_params(axis="y",direction="in")

    axs[5,3].pcolor(ϵ./ϵˢ.-1,(b./bˢ.-1)',log10.(objϵb),vmax=1.5)
    axs[5,3].set_yticklabels([])
    axs[5,3].set_xlabel("σ")
    axs[5,3].set_xticklabels([" "," ","0.0","0.1"])
    axs[5,3].tick_params(axis="x",direction="in")
    axs[5,3].tick_params(axis="y",direction="in")

    axs[5,4].pcolor(E./Eˢ.-1,(b./bˢ.-1)',log10.(objEb),vmax=1.5)
    axs[5,4].set_yticklabels([])
    axs[5,4].set_xlabel("ϵassoc")
    axs[5,4].set_xticklabels([" "," ","0.0","0.1"])
    axs[5,4].tick_params(axis="x",direction="in")
    axs[5,4].tick_params(axis="y",direction="in")

    axs[1,2].remove()
    axs[1,3].remove()
    axs[1,4].remove()
    axs[1,5].remove()
    axs[2,3].remove()
    axs[2,4].remove()
    axs[2,5].remove()
    axs[3,4].remove()
    axs[3,5].remove()
    axs[4,5].remove()


    plt.savefig(species*".pdf")
    plt.clf()
end

Nₛ = 1000
Nᵢ = 1000

θ₀ = [e.model.params.segment.values[1],e.model.params.epsilon.values[1],e.model.params.sigma.values[1]/1e-10]
σₚ = 0.016*θ₀

dₙ = Normal.(θ₀,σₚ)

θᵢ = θ₀.*ones(Nₛ,3)

dₛ = Chisq.(Nexp)

