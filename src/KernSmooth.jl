module KernSmooth

    export locpoly

    const libkernsmooth = joinpath(Pkg.dir("KernSmooth"), "deps", "libkernsmooth.so")

    include("utils.jl")
    include("locpoly.jl")

end
