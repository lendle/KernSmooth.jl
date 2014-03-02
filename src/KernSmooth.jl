module KernSmooth

    export locpoly, dpill

    const libkernsmooth = joinpath(Pkg.dir("KernSmooth"), "deps", "libkernsmooth.so")

    include("utils.jl")
    include("locpoly.jl")
    include("dpill.jl")

end
