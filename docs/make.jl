push!(LOAD_PATH,"../src/")
using Documenter,OpenSAFT

makedocs(sitename = "OpenSAFT.jl",
format = Documenter.HTML(
    # Use clean URLs, unless built as a "local" build
    canonical = "https://juliadocs.github.io/OpenSAFT.jl/dev/",
    assets = ["assets/logo.ico"],
),
    authors = "Pierre J. Walker, Hon Wa Yew and Andrés Riedemann.",
    pages = [
        "Home" => "index.md",
        "Background" => "theory/background.md",
        "User guide" => Any["Definitions"=>"user_guide/definitions.md",
                            "Basic Usage"=>"user_guide/basic_usage.md",
                            "Custom Databases"=>"user_guide/custom_dtb.md",
                            "Custom Equations of State"=>"user_guide/custom_eos.md"],
        "To-do list" => "to-do_list.md"])
