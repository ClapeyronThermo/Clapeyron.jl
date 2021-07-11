push!(LOAD_PATH,"../src/")
using Documenter,Clapeyron

makedocs(sitename = "Clapeyron.jl",
format = Documenter.HTML(
    # Use clean URLs, unless built as a "local" build
    canonical = "https://ypaul21.github.io/Clapeyron.jl/",
    assets = ["assets/logo.ico"],
    repo="github.com/ypaul21/Clapeyron.jl"
),
    authors = "Pierre J. Walker, Hon Wa Yew and Andrés Riedemann.",
    pages = [
        "Home" => "index.md",
        "Background" => "theory/background.md",
        "User guide" => Any["Basic Usage"=>"user_guide/basic_usage.md",
                            "Custom Databases"=>"user_guide/custom_dtb.md",
                            "Custom Methods"=>"user_guide/custom_methods.md",
                            "Custom Models"=>"user_guide/custom_eos.md"],
        "To-do list" => "to-do_list.md"])

        deploydocs(;
    repo="github.com/ypaul21/Clapeyron.jl.git",
    devbranch = "development",
)