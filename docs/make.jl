using AutoDiff_md
using Documenter

makedocs(;
    modules=[AutoDiff_md],
    authors="Braden Kelly",
    repo="https://github.com/BradenDKelly/AutoDiff_md.jl/blob/{commit}{path}#L{line}",
    sitename="AutoDiff_md.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://BradenDKelly.github.io/AutoDiff_md.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/BradenDKelly/AutoDiff_md.jl",
)
