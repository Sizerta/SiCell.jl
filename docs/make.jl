using Documenter
using SiCell

makedocs(
    sitename="SiCell.jl",
    remotes=nothing,
    format=Documenter.HTML(
        prettyurls=get(ENV, "CI", nothing) == "true",
        canonical="https://yourusername.github.io/SiCell.jl/stable/",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Tutorials" => "tutorials.md",
        "Visualizations" => "Visualizations.md",
        "Case Studies" => [
            "Breast Cancer Analysis" => "case_studies/case_study_breast_cancer.md",
            "PBMC3K Workflow" => "case_studies/case_study_pbmc3k.md",
            "Revealing Hidden Cellular Plasticity in Glioblastoma with TUF" => "case_studies/case_study_glioblastoma.md"
        ]
    ],
)

deploydocs(
    repo="github.com/yourusername/SiCell.jl.git",
    devbranch="main",
)