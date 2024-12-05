using Documenter
using CHMMAIRRa

makedocs(
    sitename = "CHMMAIRRa.jl",
    format = Documenter.HTML(),
    modules = [CHMMAIRRa],
    pages = [
        "Overview" => "index.md",
        "API" => "api.md",
    ]
)

deploydocs(
    repo = "github.com/MurrellGroup/CHMMAIRRa.jl.git",
    devbranch = "main",
    push_preview = true
)