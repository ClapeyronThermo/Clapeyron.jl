using Clapeyron

function density(model::EoSModel,p,T)
    return mass_density(model,p,T,[1.,1.])
end

model = SAFTgammaEMie([],[("BMIMBF4",["BMIM"=>1,"BF4"=>1])],
                         [("BMIM",["CH3"=>2,"CH2"=>3,"IM"=>1]),("BF4",["BF4"=>1])];userlocations=["IL_like.csv"],rsp_userlocations=["Salt_like.csv"])

toestimate = [
    Dict(
        :param => :epsilon,
        :indices => (3,3),
        :recombine => true,
        :lower => 100.,
        :upper => 300.,
        :guess => 500.
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
        :lower => 3.3,
        :upper => 4.5,
        :guess => 3.7
    ),
    Dict(
        :param => :sigma,
        :indices => (4,4),
        :recombine => true,
        :factor => 1e-10,
        :lower => 3.3,
        :upper => 4.5,
        :guess => 3.7
    )
]

e = Estimation(model,toestimate,["density.csv"],[:vrmodel])

optimize!(e,Clapeyron.Metaheuristics.ECA(options = Clapeyron.Metaheuristics.Options(iterations=100)))