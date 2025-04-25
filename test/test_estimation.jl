function bubble_point(model::EoSModel,T,x)
    return 1,1
end

@testset "estimation" begin
    model = SAFTgammaMie(["ethanol","water"],epsilon_mixing = :hudsen_mccoubrey)
    
    toestimate = [Dict(
        :param => :epsilon,
        :indices => (1,1),
        :recombine => true,
        :lower => 100.,
        :upper => 300.,
        :guess => 250.
    ),
    Dict(
        :param => :epsilon,
        :indices => (2,3),
        :symmetric => false,
        :lower => 100.,
        :upper => 300.,
        :guess => 250.
    ),
    Dict(
        :param => :sigma,
        :indices => (1,1),
        :factor => 1e-10,
        :lower => 3.2,
        :upper => 4.0,
        :guess => 3.7
        ),
    Dict(
        :param => :epsilon_assoc,
        :indices => 2,
        :cross_assoc => true,
        :lower => 2000.,
        :upper => 4000.,
        :guess => 3500.
        )]

    estimator,objective,initial,upper,lower = Estimation(model,toestimate,["../examples/data/bubble_point.csv"],[:vrmodel])

    model2 = return_model(estimator,model,initial)
    @test model2.params.epsilon[1,1] == initial[1] # Test that the parameter was correctly updated
    @test model2.params.epsilon[1,2] ≈ 251.57862124740765 rtol = 1e-6 # Test that the combining rule was used
    @test model2.params.epsilon[1,2] == model2.params.epsilon[1,2] # Test that the unlike parameters remain symmetric

    @test model2.params.epsilon[2,3] == initial[2] # Test that the parameter was updated
    @test model2.params.epsilon[2,3] != model2.params.epsilon[3,2] # Test that the parameter is no longer symmetric

    @test model2.params.sigma[1,1] == initial[3]*1e-10 # Test that the factor was used to update the parameters

    @test model2.params.epsilon[2,3] == initial[2] # Test that the parameter was updated
    @test model2.params.epsilon[2,3] != model2.params.epsilon[3,2] # Test that the parameter is no longer symmetric

    @test model2.params.epsilon_assoc.values.values[2] == initial[4] # Test that the association parameter was updated
    @test model2.params.epsilon_assoc.values.values[2] == model2.params.epsilon_assoc.values.values[3] # Test that the cross-association parameter was updated

    @test objective(initial) ≈ 2.4184631612655836 rtol = 1e-6

    estimator,objective,initial,upper,lower = Estimation(model,toestimate,[(2.,"../examples/data/bubble_point.csv"),(1.,"../examples/data/bubble_point.csv")],[:vrmodel])

    @test objective(initial) ≈ 7.255389483796751 rtol = 1e-6

    #error found during #365
    modelvec = Clapeyron.EoSVectorParam(model)
    @test Clapeyron.promote_model(BigFloat,modelvec) isa Clapeyron.EoSVectorParam
end
