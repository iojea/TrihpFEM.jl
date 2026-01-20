using TrihpFEM
using Documenter

DocMeta.setdocmeta!(TrihpFEM, :DocTestSetup, :(using TrihpFEM); recursive = true)

makedocs(;
    modules = [TrihpFEM],
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    authors = "Ignacio Ojea <iojea@dm.uba.ar>",
    sitename = "TrihpFEM.jl",
    # format = Documenter.HTML(;
    #     edit_link = "main",
    #     assets = String[],
    # ),
    pages = [
        "Home" => "index.md",
        "User Guide" => [
            "Mathematical background" => "basics.md",
            "Meshes" => "meshes.md",
            "Defining a problem" => "example.md"
        ],
    ],
    warnonly = [:cross_references, :missing_docs],
)
deploydocs(
    repo = "github.com/iojea/TrihpFEM.jl.git",
) 
