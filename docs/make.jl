push!(LOAD_PATH,"../src/")
using Documenter,Clapeyron,ClapeyronHANNA
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
        "Installation"=>"installation.md",
        "Tutorials" => ["Getting Started - Model Construction"=>"tutorials/basics_model_construction.md",
                        "User-defined Parameters"=>"tutorials/user_defined_parameters.md",
                        "Bulk Properties"=>"tutorials/bulk_properties.md",
                        "Mixing and Excess Functions"=>"tutorials/mixing_functions.md",
                        "Pure Saturation Properties"=>"tutorials/pure_saturation_curves.md",
                        "Binary Phase Diagrams"=>"tutorials/binary_phase_diagram.md",
                        "Ternary Phase Diagrams"=>"tutorials/ternary_phase_diagrams.md",
                        "SLE Phase Diagrams"=>"tutorials/sle_phase_diagrams.md"],
        # "How-to Guides" => ["Implement your own Methods"=>"user_guide/custom_methods.md",
        #                     "Implement your own Models"=>"user_guide/custom_model.md"],
        "Notebook Examples" => "notebook_examples.md",
        "Available EoS" => ["Ideal Models" => "eos/ideal.md"
                            "Cubic Models" => "eos/cubic.md"
                            "Activity Models" => "eos/activity.md"
                            "SAFT and CPA Models"  => "eos/saft.md"
                            "Empiric Helmholtz Models" => "eos/empiric.md"
                            "Property Correlations" =>  "eos/correlations.md"
                            "Other Models" => "eos/misc.md"],
        "Frequently Asked Questions" => "faq.md",
        "Available Properties" => ["Basic Properties" => "properties/basic.md",
                                   "Bulk Properties" => "properties/bulk.md",
                                   "Single phase Properties" => "properties/single.md",
                                   "Multiphase Properties"  => "properties/multi.md",
                                   "Electrolyte Properties" => "properties/electrolytes.md"],

        "API" => ["Parameters" => "api/parameters.md",
                  "Macros" => "api/macros.md",
                  "Association" => "api/association.md",
                  "Parameter Estimation" => "api/estimation.md"],
        "To-do list" => "to-do_list.md"]
                  )

        deploydocs(;
    repo="github.com/ClapeyronThermo/Clapeyron.jl.git",
)
