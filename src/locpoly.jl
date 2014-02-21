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

    # @assert isa(xcounts, Vector{T})
    # @assert isa(ycounts, Vector{T})
    # @assert isa(drv, Int)
    # @assert isa(delta, T)
    # @assert isa(hdisc, Vector{T})
    # @assert isa(Lvec, Vector{Int})
    # @assert isa(indic, Vector{Int})
    # @assert isa(midpts, Vector{Int})
    # @assert isa(M, Int)
    # @assert isa(Q, Int)
    # @assert isa(fkap, Vector{T})
    # @assert isa(pp, Int)
    # @assert isa(ppp, Int)
    # @assert isa(ss, Matrix{T})
    # @ascii(::Array{Uint8, 1})sert isa(tt, Matrix{T})
    # @assert isa(Smat, Matrix{T})
    # @assert isa(Tvec, Vector{T})
    # @assert isa(ipvt, Vector{Int})
    # @assert isa(curvest, Vector{T})

    ## Call FORTRAN routine "locpol"

    mid = Lvec[1] + 1
    for i in 1:Q-1
        midpts[i] = mid
        fkap[mid] = 1.0
        for j = 1:Lvec[i]
            fkap[mid+j] = exp(-(delta*j/hdisc[i])^2/2)
            fkap[mid-j] = fkap[mid+j]
        end
        mid += Lvec[i] + Lvec[i+1] + 1
    end

    midpts[Q] = mid
    fkap[mid] = 1.0
    for j = 1:Lvec[Q]
        fkap[mid+j] = exp(-(delta*j/hdisc[Q])^2/2)
        fkap[mid-j] = fkap[mid+j]
    end

    ## Combine kernel weights and grid counts
    for k = 1:M
        if xcounts[k] != 0.0
            for i in 1:Q
                for j in max(1, k-Lvec[i]):min(M, k+Lvec[i])
                    if indic[j] == i
                        fac = 1.0
                        fkap_kmjmpi = fkap[k-j+midpts[i]]
                        @inbounds ss[j, 1] += xcounts[k]*fkap_kmjmpi
                        @inbounds tt[j, 1] += ycounts[k]*fkap_kmjmpi
                        for ii = 2:ppp
                            fac *= delta * (k-j)
                            @inbounds ss[j, ii] += xcounts[k]*fkap_kmjmpi*fac
                            if ii <= pp
                                @inbounds tt[j, ii] += ycounts[k]*fkap_kmjmpi*fac
                            end
                        end
                    end
                end
            end
        end
    end

    for k in 1:M
        for i in 1:pp
            for j in 1:pp
                indss = i + j - 1

                Smat[i,j] = ss[k,indss]
            end
            Tvec[i] = tt[k,i]
        end

        # call dgefa(Smat,ipp,ipp,ipvt,info)
        # call dgesl(Smat,ipp,ipp,ipvt,Tvec,0)

        #Tvec = (Smat \ Tvec)

        curvest[k] = (Smat \ Tvec)[drv+1]
    end

    curvest = gamma(drv+1) .* curvest

    (gpoints, curvest)
end
locpoly{T<:FloatingPoint}(x::Vector{T}, bandwidth::Union(T, Vector{T}); args...) = locpoly(x, Float64[], bandwidth, args...)
