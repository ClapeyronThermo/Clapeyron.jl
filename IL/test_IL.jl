using Clapeyron

model = SAFTgammaEMie([], # Solvent
                      [("BMIMBF4",["BF4"=>1,"BMIM"=>1]),
                       ("EMIMBF4",["BF4"=>1,"EMIM"=>1])], 
                      [("BMIM",["CH3"=>2,"CH2"=>3,"IM"=>1]),
                       ("EMIM",["CH3"=>2,"CH2"=>1,"IM"=>1]),
                       ("BF4",["BF4"=>1])];
                       idealmodel=JobackIdeal,
                      userlocations=["IL_like.csv","IL_assoc.csv"], # Parameters for the groups
                      ideal_userlocations=["Ideal.csv"],
                      rsp_userlocations=["Salt_like.csv"]) # Parameters for the dielectric constant
toestimate = [
    Dict(
        :param => :epsilon,
        :indices => (3,3),
        :recombine => true,
        :lower => 100.,
        :upper => 800.,
        :guess => 200.
    ),
    Dict(
        :param => :epsilon,
        :indices => (4,4),
        :recombine => true,
        :lower => 100.,
        :upper => 800.,
        :guess => 500.
    ),
    Dict(
        :param => :sigma,
        :indices => (3,3),
        :recombine => true,
        :factor => 1e-10,
        :lower => 3.0,
        :upper => 4.5,
        :guess => 3.7
    ),
    Dict(
        :param => :sigma,
        :indices => (4,4),
        :recombine => true,
        :factor => 1e-10,
        :lower => 3.0,
        :upper => 4.5,
        :guess => 3.7
    ),
    Dict(
        :param => :shapefactor,
        :indices => (3),
        :lower => 0.5,
        :upper => 1.5,
        :guess => 1.0
    ),
    Dict(
        :param => :epsilon_assoc,
        :lower => 700.,
        :upper => 3000.,
        :guess => 2500.
    )
]

function density(model::EoSModel,p,T,id)
    z = ones(length(model.ions)).*1e-30
    z[1] = 1
    z[Int(id)+1] = 1
    return mass_density(model,p,T,z)
end

function heat_capacity(model::EoSModel,p,T,id)
    z = ones(length(model.ions)).*1e-30
    z[1] = 1
    z[Int(id)+1] = 1
    return isobaric_heat_capacity(model,p,T,z)
end

e = Estimation(model,toestimate,["density.csv"],[:vrmodel])

optimize!(e,Clapeyron.Metaheuristics.ECA(options = Clapeyron.Metaheuristics.Options(iterations=100)))

model = e.model

pred = density.(model,e.data[1].inputs[1],e.data[1].inputs[2],e.data[1].inputs[3])

println(pred)