function pH(model::ReactiveAqModel,p,T,n0;z0=nothing)
    n = equilibrate(model, p, T, n0; z0 = z0)
    a = aqueous_activity(model.eosmodel,p,T,n)
    hydronium_id = find_hydronium_index(model)
    return -log10(a[hydronium_id])
end
