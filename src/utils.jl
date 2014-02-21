## For application of linear binning to a univariate data set.
function linbin(X::Vector{Float64}, gpoints::Vector{Float64}, truncate::Bool = true)
    n = length(X)
    M = length(gpoints)

    a = gpoints[1]
    b = gpoints[M]

    trun = truncate? 1: 0

    res = Array(Float64, M)

    ccall((:linbin_, libkernsmooth), Ptr{Void},
             (Ptr{Float64}, Ptr{Int}, Ptr{Float64}, Ptr{Float64}, Ptr{Int}, Ptr{Int}, Ptr{Float64}),
             X, &n, &a, &b, &M, &trun, res)

    res
end


## For application of linear binning to a regression
## data set.
function rlbin(X::Vector{Float64}, Y::Vector{Float64}, gpoints::Vector{Float64}, truncate::Bool = true)
    n = length(X)
    M = length(gpoints)
    a = gpoints[1]
    b = gpoints[M]

    trun = truncate ? 1 : 0

    xcnts = zeros(M)
    ycnts = zeros(M)

    ccall((:rlbin_, libkernsmooth), Ptr{Void},
      (Ptr{Float64}, Ptr{Float64}, Ptr{Int}, Ptr{Float64}, Ptr{Float64}, Ptr{Int}, Ptr{Int}, Ptr{Float64}, Ptr{Float64}),
      X, Y,  &n, &a, &b, &M, &trun, xcnts, ycnts)

    (xcnts, ycnts)
end
