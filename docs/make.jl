push!(LOAD_PATH,"../src/")
using Documenter,Clapeyron

makedocs(sitename = "Clapeyron.jl",
format = Documenter.HTML(
    # Use clean URLs, unless built as a "local" build
    canonical = "https://ClapeyronThermo.github.io/Clapeyron.jl/",
    assets = ["assets/logo.ico"],
),
warnonly = Documenter.except(),
    authors = "Pierre J. Walker, Hon Wa Yew and Andrés Riedemann.",
    pages = [
        "Home" => "index.md",
        "Background" => "theory/background.md",
        "Basic Usage" => "user_guide/basic_usage.md",
        "Customization" => Any[
                            "Custom Databases"=>"user_guide/custom_dtb.md",
                            "Custom Methods"=>"user_guide/custom_methods.md",
                            "Custom Models"=>"user_guide/custom_model.md"],
        "Notebook Examples" => "notebook_examples.md",
        "To-do list" => "to-do_list.md",
        
        "Available EoS" => [
        "Ideal Models" => "eos/ideal.md"
        "Cubic Models" => "eos/cubic.md"
        "Activity Models" => "eos/activity.md"
        "SAFT and CPA Models"  => "eos/saft.md"
        "Electrolyte Models" => "eos/electrolytes.md"
        "Empiric Helmholtz Models" => "eos/empiric.md"
        "Property Correlations" =>  "eos/correlations.md"
        "Other Models" => "eos/misc.md"
        ],
        
        "Available Properties" => [
            "Basic Properties" => "properties/basic.md",
            "Bulk Properties" => "properties/bulk.md",
            "Single phase Properties" => "properties/single.md",
            "Multiphase Properties"  => "properties/multi.md",
            "Electrolyte Properties" => "properties/electrolytes.md",
            ],

        "API" => Any[
        "Parameters" => "api/parameters.md",
        "Macros" => "api/macros.md",
        "Association" => "api/association.md",
        "Parameter Estimation" => "api/estimation.md",
        ],
        "Developer Guide" => "dev.md"])

        deploydocs(;
    repo="github.com/ClapeyronThermo/Clapeyron.jl.git",
)
