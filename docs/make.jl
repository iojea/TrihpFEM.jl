using TrihpFEM
using Documenter

DocMeta.setdocmeta!(TrihpFEM, :DocTestSetup, :(using TrihpFEM); recursive=true)

makedocs(;
    modules=[TrihpFEM],
    authors="Ignacio Ojea <iojea@dm.uba.ar>",
    sitename="TrihpFEM.jl",
    format=Documenter.HTML(;
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
    warnonly = [:cross_references, :missing_docs],
)
