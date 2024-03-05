using MyFirstPackage
using Documenter

DocMeta.setdocmeta!(MyFirstPackage, :DocTestSetup, :(using MyFirstPackage); recursive=true)

makedocs(;
    modules=[MyFirstPackage],
    authors="zhangdezheng",
    sitename="MyFirstPackage.jl",
    format=Documenter.HTML(;
        canonical="https://zhangdezheng.github.io/MyFirstPackage.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/zhangdezheng/MyFirstPackage.jl",
    devbranch="main",
)
