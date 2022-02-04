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
        "User guide" => Any["Basic Usage"=>"user_guide/basic_usage.md",
                            "Custom Databases"=>"user_guide/custom_dtb.md",
                            "Custom Methods"=>"user_guide/custom_methods.md",
                            "Custom Models"=>"user_guide/custom_model.md"],
        "To-do list" => "to-do_list.md",
        
        "Available EoS" => [
        "Ideal Models" => "eos/ideal.md"
        "Cubic Models" => "eos/cubic.md"
        "Empiric Helmholtz Models" => "eos/empiric.md"
        ],
        
        "API" => Any[
        "Parameters" => "api/parameters.md",
        "Macros" => "api/macros.md",
        "Properties" => "api/properties.md",
        "Automatic Differenciation" => "api/ad.md",
        ]
        ])

        deploydocs(;
    repo="github.com/ypaul21/Clapeyron.jl.git",
)