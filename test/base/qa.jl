#File to add any quality assurance tests.

@printline

@testset verbose = true "QA" begin
    test_ambiguities(Clapeyron)
    test_unbound_args(Clapeyron)
    test_undefined_exports(Clapeyron)
    test_project_extras(Clapeyron)
    #test_stale_deps(Clapeyron)
    #test_deps_compat(Clapeyron)
    test_piracies(Clapeyron)
    test_persistent_tasks(Clapeyron)
    test_undocumented_names(Clapeyron)
end
