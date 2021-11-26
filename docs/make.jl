using Documenter
using DDA


makedocs(
    sitename = "DDA",
    authors="Adam Fekete <adam.fekete@unamur.be> and contributors",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical="https://fekad.github.io/DDA.jl"
    ),
    modules = [DDA],
    pages=[
        "Home" => "index.md",
        "Background" => "theory.md",
        "Polarizability" => "polarizability.md"
    ]
)


deploydocs(
    repo = "github.com/fekad/DDA.jl",
    devbranch = "main"
)