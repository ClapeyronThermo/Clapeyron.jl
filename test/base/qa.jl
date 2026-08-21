#File to add any quality assurance tests.

@printline

@testset verbose = true "QA" begin
    Aqua.test_ambiguities(Clapeyron,broken = true)
    Aqua.test_unbound_args(Clapeyron,broken = true)
    Aqua.test_undefined_exports(Clapeyron)
    Aqua.test_project_extras(Clapeyron)
    #test_stale_deps(Clapeyron)
    #test_deps_compat(Clapeyron)
    Aqua.test_piracies(Clapeyron)
    Aqua.test_persistent_tasks(Clapeyron)
    Aqua.test_undocumented_names(Clapeyron,broken = true)
end
