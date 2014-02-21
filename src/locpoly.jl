## For computing a binned local polynomial
## regression estimator of a univariate regression
## function or its derivative.
## The data are discretised on an equally
## spaced grid. The bandwidths are discretised on a
## logarithmically spaced grid.
function locpoly{T<:FloatingPoint}(x::Vector{T}, y::Vector{T}, bandwidth::Union(T, Vector{T});
    drv::Int = 0,
    degree::Int=drv+1,
    kernel::Symbol = :normal,
    gridsize::Int = 401,
    bwdisc::Int = 25,
    range_x::Vector{T}=T[],
    binned::Bool = false,
    truncate::Bool = true)

    if isa(bandwidth, T) && bandwidth <  0.0 
        error("'bandwidth' must be strictly positive")
    end

    if range_x == T[] && !binned
        maxx = maximum(x)
        minx = minimum(x)
        if y == T[]
            extra = 0.05 * (maxx - minx)
            range_x = [minx - extra, maxx + extra]
        else
            range_x = [minx, maxx]
        end
    end

    ## Rename common variables
    M = gridsize
    Q = bwdisc
    a = range_x[1]
    b = range_x[2]
    pp = degree + 1
    ppp = 2*degree + 1
    tau = 4

    ## Decide whether a density estimate or regression estimate is required.

    if y == T[]    # obtain density estimate
        n = length(x)
        gpoints = linspace(a, b, M)
        xcounts = linbin(x, gpoints, truncate)
        ycounts = (M-1) .* xcounts ./ (n*(b-a))
        xcounts = ones(M) #rep(1, M)
    else           # obtain regression estimate
        ## Bin the data if not already binned
        if !binned
            gpoints = linspace(a, b, M)
            xcounts, ycounts = rlbin(x, y, gpoints, truncate)
        else
            xcounts = x
            ycounts = y
            M = length(xcounts)
            gpoints = linspace(a, b, M)
        end
    end

    ## Set the bin width
    delta = (b-a)/(M-1)

    indc = T[]
    Lvec = Int[]
    hdisc = T[]
    indic = Int[]
    ## Discretise the bandwidths
    if length(bandwidth) == M
        hlow = minimum(bandwidth)
        hupp = maximum(bandwidth)
        hdisc = [exp(h) for h in linspace(log(hlow),log(hupp),Q)]

        ## Determine value of L for each member of "hdisc"
        Lvec = [ifloor(tau*h/delta) for h in hdisc]

        ## Determine index of closest entry of "hdisc"
        ## to each member of "bandwidth"
        indic = if Q > 1
                    gap = (log(hdisc[Q])-log(hdisc[1]))/(Q-1)
                    if gap == 0
                        ones(Int, M)
                    else
                        lhlow = log(hlow)
                        [iround((log(b) - lhlow)/gap + 1.0) for b in bandwidth]
                    end
                else
                    ones(Int, M)
                end
    elseif length(bandwidth) == 1
        indic = ones(Int, M)
        Q = 1
        Lvec = fill(ifloor(tau*bandwidth/delta), Q)
        hdisc = fill(bandwidth, Q)
    else
        error("'bandwidth' must be a vector of length one or of length 'gridsize'")
    end

    if minimum(Lvec) == 0
        error("Binning grid too coarse for current (small) bandwidth: consider increasing 'gridsize'")
    end

    ## Allocate space for the kernel vector and final estimate

    dimfkap = 2 * sum(Lvec) + Q
    fkap = zeros(dimfkap)
    curvest = zeros(M)
    midpts = zeros(Int, Q)
    ss = zeros(M, ppp)
    tt = zeros(M, pp)
    Smat = zeros(pp, pp)
    Tvec = zeros(pp)
    ipvt = zeros(Int, pp)


    ccall((:locpol_, libkernsmooth), Ptr{Void},
      (Ptr{Float64}, Ptr{Float64}, Ptr{Int}, Ptr{Float64}, Ptr{Float64}, Ptr{Int}, Ptr{Int}, Ptr{Int},
       Ptr{Int}, Ptr{Int}, Ptr{Float64}, Ptr{Int}, Ptr{Int}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
       Ptr{Float64}, Ptr{Int}, Ptr{Float64}),
      xcounts, ycounts,  &drv, &delta, hdisc, Lvec, indic, midpts, &M, &Q, fkap, &pp, &ppp, 
      ss, tt, Smat, Tvec, ipvt, curvest)


    curvest = gamma(drv+1) .* curvest

    (gpoints, curvest)
end

locpoly{T<:FloatingPoint}(x::Vector{T}, bandwidth::Union(T, Vector{T}); args...) = locpoly(x, Float64[], bandwidth, args...)
