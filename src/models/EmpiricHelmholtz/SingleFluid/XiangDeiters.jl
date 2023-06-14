function XiangDeitersConsts(ω,θ)
    d = [1, 1, 1, 2, 3, 7, 1, 1, 2, 5, 1, 1, 4, 2]
    t = [0.25, 1.25, 1.5, 1.375, 0.25, 0.875, 0, 2.375, 2, 2.125, 3.5, 6.5, 4.75, 12.5]
    l = [1, 1, 1, 1, 2, 2, 2, 3]
    a0 = [8.5740489E-01,  -3.2863233E+00, 1.6480939E+00,  -5.4524817E-02, 6.1623592E-02, 2.7389266E-04,  -6.0655087E-02,
    -3.1811852E-02, -1.1550422E-01, -1.8610466E-02, -1.8348671E-01, 5.5071325E-03, -1.2268039E-02, -5.0433436E-03]
    a1 = [5.6200117E-01, 3.2439544E+00, -4.9628768E+00, -2.2132851E-01, 9.3388356E-02, 2.4940171E-05,  -1.7519204E-01,
    8.9325660E-01, 2.9886613E+00, 1.0881387E-01,  -6.7166746E-01, 1.4477326E-01, -2.8716809E-01, -1.1478402E-01]
    a2 = [-8.1680511E+01, 4.6384732E+02,  -2.7970850E+02, 2.9317364E+01, -2.2324825E+01, -5.0932691E-02, -7.2836590E+00,
    -2.2063100E+02, -3.0435126E+02, 5.8514719E+00,  1.7995451E+02, -1.0178400E+02, 4.0848053E+01,  1.2411984E+01]
    n = a0 .+ a1 .* ω + a2 .* θ
    return EmpiricSingleFluidResidualParam(n,t,d,l)
end

function XiangDeiters(components;
    idealmodel = BasicIdeal,
    userlocations = String[],
    ideal_userlocations = String[],
    Rgas = nothing,
    verbose = false)
    params = getparams(components, ["properties/critical.csv", "properties/molarmass.csv"]; userlocations=userlocations, verbose=verbose)
    Pc = params["Pc"][1]
    Vc = params["Vc"][1]
    Mw = params["Mw"][1]
    Tc = params["Tc"][1]
    acentricfactor = params["acentricfactor"][1]

    Zc = Pc*Vc/(R̄*Tc)
    θ = (Zc - 0.29)^2
    rhoc = 1/Vc
    lb_volume = 1/(3.25*rhoc)
    residual = XiangDeitersConsts(acentricfactor,θ)
    if Rgas === nothing
        R = R̄
    else
        R = Rgas
    end
    properties = ESFProperties(Mw,Tc,rhoc,lb_volume,Tc,Pc,rhoc,NaN,NaN,NaN,NaN,acentricfactor,Rgas)
    ancillaries = propane_ancillary_cs(components,Tc,Pc,Vc)
    init_idealmodel = init_model(idealmodel,components,ideal_userlocations,false)
    ideal = EmpiricSingleFluidIdealParamBuilder(init_idealmodel,properties)
    EmpiricSingleFluid(components,properties,ancillaries,ideal,residual,String["10.1016/j.ces.2007.11.029"])
end

function EmpiricSingleFluidIdealParamBuilder(model::EoSModel,props)
    return _parse_ideal(ideal_dict_format(model,props))    
end

function ideal_dict_format(model::BasicIdealModel,props)
    [
        Dict(:type => "IdealGasHelmholtzLead",
        :a1 => -1.0,
        :a2 => 0.0
        ),

        Dict(:type => "IdealGasHelmholtzLogTau",
        :a => -1.5
        ),
    ]
end

function ideal_dict_format(model::ReidIdealModel,props)
    Tc = props.Tc
    [
        Dict(:type => "IdealGasHelmholtzCP0PolyT",
        :c => [model.params.coeffs[1]...]
        :t => [0,1,2,3]
        :Tc => Tc
        :T0 => 298.15
        )
    ]
end
