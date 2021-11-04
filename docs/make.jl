using PCCA
using Documenter

DocMeta.setdocmeta!(PCCA, :DocTestSetup, :(using PCCA); recursive=true)

makedocs(;
    modules=[PCCA],
    authors="Alexander Sikorski",
    repo="https://github.com/axsk/PCCA.jl/blob/{commit}{path}#{line}",
    sitename="PCCA.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://axsk.github.io/PCCA.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/axsk/PCCA.jl",
    devbranch="main",
)
