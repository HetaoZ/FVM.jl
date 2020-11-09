using FVM
using Documenter

makedocs(;
    modules=[FVM],
    authors="Hetao Z.",
    repo="https://github.com/HetaoZ/FVM.jl/blob/{commit}{path}#L{line}",
    sitename="FVM.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://HetaoZ.github.io/FVM.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/HetaoZ/FVM.jl",
)
