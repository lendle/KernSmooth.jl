## For application of linear binning to a univariate data set.
function linbin{T<:FloatingPoint}(X::Vector{T}, gpoints::Vector{T}, truncate::Bool = true)
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
function rlbin{T<:FloatingPoint}(X::Vector{T}, Y::Vector{T}, gpoints::Vector{T}, truncate::Bool = true)
    n = length(X)
    M = length(gpoints)
    a = gpoints[1]
    b = gpoints[M]

    xcnts = zeros(M)
    ycnts = zeros(M)

    delta = (b-a)/(M - 1.0)

    for i in 1:length(X)
        lxi = ((X[i]-a)/delta) + 1

        li = ifloor(lxi)
        rem = lxi - li
        if  1 <= li < M
            xcnts[li] += 1-rem
            xcnts[li+1] += rem
            ycnts[li] += (1-rem) * Y[i]
            ycnts[li+1] += rem * Y[i]
        end

        if (li < 1 && !truncate)
            xcnts[1] += 1.0
            ycnts[1] += Y[i]
        end

        if li >= M && !truncate
            xcnts[M] += 1.0
            ycnts[M] += Y[i]
        end
    end

    (xcnts, ycnts)
end
