using Pkg; Pkg.activate("./docs")

using Documenter, GCATConductance
makedocs(modules = [GCATConductance], doctest = true, sitename = "GCAT-Conductance", remotes = nothing)