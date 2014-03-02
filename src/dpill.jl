## Computes a direct plug-in selector of the
## bandwidth for local linear regression as
## described in the 1996 J. Amer. Statist. Assoc.
## paper by Ruppert, Sheather and Wand.

function dpill(x::Vector{Float64}, y::Vector{Float64};
               blockmax::Int = 5,
               divisor::Int = 20,
               trim::Float64 = 0.01,
               proptrun::Float64 = 0.05,
               gridsize::Int = 401,
               range_x::Vector{Float64} = Float64[],
               truncate = true)
    ## Trim the 100(trim)% of the data from each end (in the x-direction).
    sp = sortperm(x)
    x = x[sp] #don't sort in place so original x and y are not muted
    y = y[sp]
    indlow = ifloor(trim*length(x)) + 1
    indupp = length(x) - ifloor(trim*length(x))

    x = x[indlow:indupp]
    y = y[indlow:indupp]

    ## Rename common parameters
    n = length(x)
    M = gridsize
    if range_x == Float64[]
        range_x = [minimum(x), maximum(x)]
    end

    a = range_x[1]
    b = range_x[2]

    ## Bin the data

    gpoints = linspace(a, b, M)
    xcounts, ycounts = rlbin(x, y, gpoints, truncate)

    ## Choose the value of N using Mallow's C_p
    Nmax = max(min(ifloor(n/divisor), blockmax), 1)
    Nval = cpblock(x, y, Nmax, 4)

    ## Estimate sig^2, theta_22 and theta_24 using quartic fits
    ## on "Nval" blocks.

    (sigsqQ, th22Q, th24Q) = blkest(x, y, Nval, 4)

    ## Estimate theta_22 using a local cubic fit
    ## with a "rule-of-thumb" bandwidth: "gamseh"

    gamseh = (sigsqQ*(b-a)/(abs(th24Q)*n))

    if th24Q < 0.0
        gamseh = (3.0 * gamseh / (8.0 * sqrt(pi)))^(1.0/7.0)
    end
    if th24Q > 0.0
        gamseh = (15.0 * gamseh / (16.0 * sqrt(pi)))^(1.0/7.0)
    end

    mddest = locpoly(xcounts, ycounts, gamseh, drv=2, range_x=range_x, binned=true)[2]

    llow = ifloor(proptrun*M) + 1
    lupp = M - ifloor(proptrun*M)
    th22kn = sum((mddest[llow:lupp].^2) .* xcounts[llow:lupp])/n

    ## Estimate sigma^2 using a local linear fit
    ## with a "direct plug-in" bandwidth: "lamseh"
    C3K = (1.0/2.0) + 2.0*sqrt(2.0) - (4.0/3.0)*sqrt(3.0)
    C3K = (4.0*C3K/(sqrt(2.0*pi)))^(1.0/9.0)
    lamseh = C3K*(((sigsqQ^2)*(b-a)/((th22kn*n)^2))^(1.0/9.0))

    ## Now compute a local linear kernel estimate of
    ## the variance.
    mest = locpoly(xcounts, ycounts, lamseh, range_x=range_x, binned=true)[2]
    Sdg = sdiag(xcounts, lamseh, range_x=range_x, binned=true)[2]
    SSTdg = sstdiag(xcounts, lamseh, range_x=range_x, binned=true)[2]
    sigsqn = sum(y.^2) .- 2*sum(mest.*ycounts) + sum((mest.^2).*xcounts)
    sigsqd = n - 2*sum(Sdg.*xcounts) + sum(SSTdg.*xcounts)
    sigsqkn = sigsqn/sigsqd

    ## Combine to obtain final answer.
    (sigsqkn*(b-a)/(2.0*sqrt(pi)*th22kn*n))^(1.0/5.0)
end
