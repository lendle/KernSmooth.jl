## For application of linear binning to a univariate data set.
function linbin{T<:FloatingPoint}(X::Vector{T}, gpoints::Vector{T}, truncate::Bool = true)
    M = length(gpoints)

    a = gpoints[1]
    b = gpoints[M]

    gcnts = zeros(M)
    delta = (b-a)/(M-1.0)
    for i in 1:length(X)
        lxi = (X[i]-a)/delta + 1

        li = ifloor(lxi)

        rem = lxi - li

        if 1 <= li < M
            gcnts[li] += (1-rem)
            gcnts[li+1] += rem
        end

        if li < 1 && !truncate
            gcnts[1] += 1.0
        end

        if li >= M && !truncate
            gcnts[M] += 1.0
        end
    end
    gcnts
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
            ycnts[li] += (1-rem) * y[i]
            ycnts[li+1] += rem * y[i]
        end

        if (li < 1 && !truncate)
            xcnts[1] += 1.0
            ycnts[1] += y[i]
        end

        if li >= M && !truncate
            xcnts[M] += 1.0
            ycnts[M] += y[i]
        end
    end

    (xcnts, ycnts)
end
