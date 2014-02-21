using KernSmooth, RDatasets, Base.Test

expected = readcsv(joinpath(Pkg.dir("KernSmooth"), "test", "geyserfit.csv"), Float64, has_header=true)[1]

dat = dataset("MASS", "geyser")

x = array(dat[:Duration])
y = convert(Vector{Float64}, array(dat[:Waiting]))


resx, resy= locpoly(x, y, 0.25)

@test_approx_eq_eps(maximum(abs(resx - expected[:,1])), 0.0, sqrt(eps(Float64)))
@test_approx_eq_eps(maximum(abs(resy - expected[:,2])), 0.0, sqrt(eps(Float64)))
