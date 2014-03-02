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

## Chooses the number of blocks for the preliminary
## step of a plug-in rule using Mallows' C_p.
function cpblock(X::Vector{Float64}, Y::Vector{Float64}, Nmax::Int, q::Int)
    n = length(X)

    ## Sort the (X, Y) data with respect to the X's.

    sp = sortperm(X)
    X = X[sp] #don't sort in place so original x and y are not muted
    Y = Y[sp]

    ## Set up arrays for FORTRAN subroutine "cp"

    qq = q + 1
    RSS = zeros(Nmax)
    Xj = zeros(n)
    Yj = zeros(n)
    coef = zeros(qq)
    Xmat = zeros(n, qq)
    Cpvals = zeros(Nmax)
    wk = zeros(n)
    qraux = zeros(qq)


    ccall((:cp_, libkernsmooth), Ptr{Void},
        (Ptr{Float64}, Ptr{Float64}, Ptr{Int}, Ptr{Int}, Ptr{Int}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
         Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}),
        X, Y, &n, &qq, &Nmax, RSS, Xj, Yj, coef, Xmat, wk, qraux, Cpvals)

    sortperm(Cpvals)[1]
end

## For obtaining preliminary estimates of
## quantities required for the "direct plug-in"
## regression bandwidth selector based on
## blocked qth degree polynomial fits.

function blkest(x::Vector{Float64}, y::Vector{Float64}, Nval::Int, q::Int)
    n = length(x)

    ## Sort the (x, y) data with respect to
    ## the x's.

    sp = sortperm(x)
    x = x[sp] #don't sort in place so original x and y are not muted
    y = y[sp]

    ## Set up arrays for FORTRAN program "blkest"

    qq = q + 1
    xj = zeros(n)
    yj = zeros(n)
    coef = zeros(qq)
    Xmat = zeros(n, qq)
    wk = zeros(n)
    qraux = zeros(qq)
    sigsqe = zeros(1)
    th22e = zeros(1)
    th24e = zeros(1)

    ccall((:blkest_, libkernsmooth), Ptr{Void},
          (Ptr{Float64}, Ptr{Float64}, Ptr{Int}, Ptr{Int}, Ptr{Int}, Ptr{Int}, Ptr{Float64},
           Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
           Ptr{Float64}, Ptr{Float64}),
          x, y, &n, &q, &qq, &Nval, xj, yj, coef, Xmat, wk, qraux, sigsqe, th22e, th24e)


    (sigsqe[1], th22e[1], th24e[1])
end


function discretise_the_bandwidths(bandwidth, M, Q, tau, delta)
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

    (indic, Q, Lvec, hdisc)
end

function sdiag(x::Vector{Float64}, bandwidth;
                drv::Int=0,
                degree::Int=1,
                kernel::Symbol = :normal,
                gridsize::Int = 401,
                bwdisc::Int = 25,
                range_x::Vector{Float64}=Float64[],
                binned::Bool = false,
                truncate::Bool = true)

    if range_x == Float64[] && !binned
        maxx = maximum(x)
        minx = minimum(x)
        if y == Float64[]
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

    ## Bin the data if not already binned

    if !binned
        gpoints = linspace(a, b, M)
        xcounts = linbin(x, gpoints, truncate)
    else
        xcounts = x
        M = length(xcounts)
        gpoints = linspace(a, b, M)
    end

    ## Set the bin width

    delta = (b-a)/(M-1.0)

    ## Discretise the bandwidths

    (indic, Q, Lvec, hdisc) = discretise_the_bandwidths(bandwidth, M, Q, tau, delta)


    dimfkap = 2 * sum(Lvec) + Q
    fkap = zeros(dimfkap)
    midpts = zeros(Int, Q)
    ss = zeros(M, ppp)
    Smat = zeros(pp, pp)
    work = zeros(pp)
    det = zeros(2)
    ipvt = zeros(Int, pp)
    Sdg = zeros(M)

    ccall((:sdiag_, libkernsmooth), Ptr{Void},
          (Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Int}, Ptr{Int}, Ptr{Int}, Ptr{Int},
           Ptr{Int}, Ptr{Float64}, Ptr{Int}, Ptr{Int}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
           Ptr{Float64}, Ptr{Int}, Ptr{Float64}),
          xcounts, &delta, hdisc, Lvec, indic, midpts, &M, &Q, fkap, &pp, &ppp, ss, Smat, work,
          det, ipvt, Sdg)

    (gpoints, Sdg)
end


## For computing the binned diagonal entries of SS^T
## where S is a smoother matrix for local polynomial
## kernel regression.

function sstdiag(x::Vector{Float64}, bandwidth;
                drv::Int=0,
                degree::Int=1,
                kernel::Symbol = :normal,
                gridsize::Int = 401,
                bwdisc::Int = 25,
                range_x::Vector{Float64}=Float64[],
                binned::Bool = false,
                truncate::Bool = true)
    if range_x == Float64[] && !binned
        maxx = maximum(x)
        minx = minimum(x)
        if y == Float64[]
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

    ## Bin the data if not already binned
    if !binned
        gpoints = linspace(a, b, M)
        xcounts = linbin(x, gpoints, truncate)
    else
        xcounts = x
        M = length(xcounts)
        gpoints = linspace(a, b, M)
    end

    ## Set the bin width

    delta = (b-a)/(M-1.0)

    ## Discretise the bandwidths
    (indic, Q, Lvec, hdisc) = discretise_the_bandwidths(bandwidth, M, Q, tau, delta)

    dimfkap = 2 * sum(Lvec) + Q
    fkap = zeros(dimfkap)
    midpts = zeros(Int, Q)
    ss = zeros(M, ppp)
    uu = zeros(M, ppp)
    Smat = zeros(pp, pp)
    Umat = zeros(pp, pp)
    work = zeros(pp)
    det = zeros(2)
    ipvt = zeros(Int, pp)
    SSTd = zeros(M)

    ccall((:sstdg_, libkernsmooth), Ptr{Void},
          (Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Int}, Ptr{Int}, Ptr{Int}, Ptr{Int},
           Ptr{Int}, Ptr{Float64}, Ptr{Int}, Ptr{Int}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
           Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Int}, Ptr{Float64}),
          xcounts, &delta, hdisc, Lvec, indic, midpts, &M, &Q, fkap, &pp, &ppp, ss, uu, Smat,
          Umat, work, det, ipvt, SSTd)

    (gpoints, SSTd)
end

