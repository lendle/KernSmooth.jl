# KernSmooth

[![Build Status](https://travis-ci.org/lendle/KernSmooth.jl.png)](https://travis-ci.org/lendle/KernSmooth.jl)

## About

KernSmooth.jl is the start of a direct port of the R package [KernSmooth](http://cran.r-project.org/web/packages/KernSmooth/index.html), (v2.23-10.)
The R package carries an unlimited license.

Currently the `locpoly` function is implemented, which uses local polynomials to estimate pdf of a single variable or a regression function for two variables, or their derivatives.

## Usage

The method signatures:
```julia
locpoly{T<:FloatingPoint}(x::Vector{T}, y::Vector{T}, bandwidth::Union(T, Vector{T});
    drv::Int = 0,
    degree::Int=drv+1,
    kernel::Symbol = :normal,
    gridsize::Int = 401,
    bwdisc::Int = 25,
    range_x::Vector{T}=T[],
    binned::Bool = false,
    truncate::Bool = true)

locpoly{T<:FloatingPoint}(x::Vector{T}, bandwidth::Union(T, Vector{T}); args...)
```

* `x` - vector of x data
* `y` - vector of y data. For density estimation (of `x`), `y` should be omitted or be an empty `Vector{T}`
* `bandwidth` - should be a scalar or vector of length `gridsize`
* Other arguments are optional. For their descriptions, see the [R documentation](https://stat.ethz.ch/R-manual/R-devel/library/KernSmooth/html/locpoly.html)

