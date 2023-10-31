using PCCAPlus
using Documenter

DocMeta.setdocmeta!(PCCAPlus, :DocTestSetup, :(using PCCAPlus); recursive=true)

makedocs(;
    modules=[PCCAPlus],
    authors="Alexander Sikorski",
    repo="https://github.com/axsk/PCCAPlus.jl/blob/{commit}{path}#{line}",
    sitename="PCCAPlus.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://axsk.github.io/PCCAPlus.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/axsk/PCCAPlus.jl",
    devbranch="main",
)
