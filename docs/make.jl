using Documenter
push!(LOAD_PATH,"../src/")
makedocs(sitename = "OpenSAFT.jl",
format = Documenter.HTML(
    # Use clean URLs, unless built as a "local" build
    canonical = "https://juliadocs.github.io/OpenSAFT.jl/dev/",
    assets = ["assets/logo.ico"],
),
    authors = "Pierre J. Walker and Hon Wa Yew.",
    pages = [
        "Home" => "index.md",
        "User guide" => Any["System"=>"user_guide/system.md",
                            "Bulk properties"=>"user_guide/bulk_props.md",
                            "VLE properties"=>"user_guide/vle_props.md"],
        "To-do list" => "to-do_list.md",
        "Theory"=>Any["Methods"=>"theory/methods.md",
                      "Models"=>"theory/models.md"]])
