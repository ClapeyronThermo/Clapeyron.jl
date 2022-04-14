push!(LOAD_PATH,"../src/")
using Documenter,Clapeyron

makedocs(sitename = "Clapeyron.jl",
format = Documenter.HTML(
    # Use clean URLs, unless built as a "local" build
    canonical = "https://ypaul21.github.io/Clapeyron.jl/",
    assets = ["assets/logo.ico"],
),
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
        "Empiric Helmholtz Models" => "eos/empiric.md"
        "Other Models" => "eos/misc.md"
        ],
        
        "API" => Any[
        "Parameters" => "api/parameters.md",
        "Macros" => "api/macros.md",
        "Properties" => "api/properties.md",
        ]
        ])

        deploydocs(;
    repo="github.com/ypaul21/Clapeyron.jl.git",
)