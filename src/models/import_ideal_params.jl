function create_IdealParams(components, raw_params, ideal_model)
  if ideal_model == "Monomer"
    ideal_params_dict = filterparams_ideal(raw_params, ["Mw"])
    Mw = ideal_params_dict["Mw"]
    map!(x->x*1E-3, values(Mw))
    return Monomer(components,MonomerParams(Mw))
  elseif ideal_model == "Walker"
    ideal_params_dict = filterparams_ideal(raw_params, ["Mw","Nrot","theta_1","theta_2","theta_3","theta_4",
                                                  "deg_1","deg_2","deg_3","deg_4"])
    Mw = ideal_params_dict["Mw"]
    Nrot = ideal_params_dict["Nrot"]

    theta_V = Dict()
    deg_V   = Dict()
    for i in components
      theta_V[i] = [ideal_params_dict["theta_1"][i],
                    ideal_params_dict["theta_2"][i],
                    ideal_params_dict["theta_3"][i],
                    ideal_params_dict["theta_4"][i]]
      deg_V[i]   = [ideal_params_dict["deg_1"][i],
                    ideal_params_dict["deg_2"][i],
                    ideal_params_dict["deg_3"][i],
                    ideal_params_dict["deg_4"][i]]
    end
    return Walker(components,WalkerParams(Mw,Nrot,theta_V,deg_V))
  elseif ideal_model == "Reid"
    ideal_params_dict = filterparams_ideal(raw_params, ["a","b","c","d"])
    poly_coef = Dict()
    for i in components
      poly_coef[i] = [ideal_params_dict["a"][i],
                      ideal_params_dict["b"][i],
                      ideal_params_dict["c"][i],
                      ideal_params_dict["d"][i]]
    end
    return Reid(components,ReidParams(poly_coef))
  elseif  ideal_model == "Basic"
    return Basic(components)
  end
end
