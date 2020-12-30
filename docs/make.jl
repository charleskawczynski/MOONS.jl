using Documenter, MOONS

mathengine = MathJax(Dict(
    :TeX => Dict(
        :equationNumbers => Dict(:autoNumber => "AMS"),
        :Macros => Dict(),
    ),
))

format = Documenter.HTML(
    prettyurls = get(ENV, "CI", "") == "true",
    mathengine = mathengine,
    collapselevel = 1,
)

makedocs(
    sitename = "MOONS.jl",
    format = format,
    doctest = false,
    strict = true,
    checkdocs = :exports,
    clean = true,
    modules = [MOONS],
    pages = Any[
        "Home" => "index.md",
        "Grids" => "grids.md",
        "API" => "api.md",
    ],
)

deploydocs(
    repo = "github.com/charleskawczynski/MOONS.jl.git",
    target = "build",
    push_preview = true,
    forcepush = true,
)
