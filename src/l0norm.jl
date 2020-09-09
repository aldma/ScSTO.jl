"""
Nonnegativity-constrained proximal operator of nnz()
    nnz_y = proxl0nonneg!( x, a, y )
solves the following problem
    min_y    a nnz(y) + 1/2 ||y-x||^2    s.t.    y ≥ 0
where `a` is a given nonnegative scalar.
The analytical expression is available:
    y = { x  if  x ≥ √2a,  0  otherwise
"""
function proxl0nonneg!( x::Vector{R}, a::R, y::Vector{R} ) where {R <: Real}
    @assert a >= R(0)
    y .= x
    if a > R(0)
        y[.!(x .> sqrt( 2 * a ))] .= R(0)
    else
        y .= max.( R(0), x )
    end
    nnz_y = sum(y .> R(0))
    return nnz_y
end

"""
Simplex-constrained proximal operator of nnz()
    nnz_y = proxl0simplex!( x, a, b, y )
solves the following problem
    min_y    a nnz(y) + 1/2 ||y-x||^2    s.t.    y ≥ 0,    1 ⋅ y = b
where `a` and `b` are given nonnegative scalars.
"""
function proxl0simplex!( x::Vector{R}, a::R, b::R, y::Vector{R} ) where {R <: Real}
    @assert a >= R(0)
    @assert b >= R(0)
    if a == R(0)
        projsimplex!(x, b, y)
        nnz_y = sum(y .> R(0))
        return nnz_y
    end
    n = length(x)
    y .= x # max.( x, 0.0 ) # project
    prm = sortperm(y) # sort
    y .= y[prm]
    mv = Vector(0:n-1) # number of zeros
    lv = (b .- reverse(cumsum(reverse(y)))) ./ (n .- mv)
    sv = cumsum([R(0); y[1:n-1]] .^ 2)
    feas = (y .+ lv .> R(0)) # feasibility
    @assert any(feas)
    mv = mv[feas]
    lv = lv[feas]
    sv = sv[feas]
    cv = a .* (n .- mv) .+ 0.5 .* sv .+ 0.5 .* (n .- mv) .* (lv .^ 2) # cost
    c, i = findmin(cv) # find minimum
    m = mv[i]
    l = lv[i]
    y[1:m] .= R(0)
    y[m+1:n] .+= l
    nnz_y = n - m
    y .= y[invperm(prm)] # unsort
    return nnz_y
end

"""
Projection onto the simplex
    projsimplex!(x, b, y)
corresponds to solving the following problem
        min_y    1/2 ||y-x||^2    s.t.    y ≥ 0,    1 ⋅ y = b
where `b` is a given nonnegative scalar.
See arxiv.org/abs/1101.6081
"""
function projsimplex!(x::Vector{R}, b::R, y::Vector{R}) where {R <: Real}
    @assert b >= R(0)
    if b == R(0)
        y .= R(0)
        return nothing
    end
    n = length(x)
    y .= sort(x, rev=true) # sort
    flag = false
    tmpsum = R(0)
    for i in 1:n-1
        tmpsum += y[i] # cumulative sum
        tmpmax = (tmpsum - b) / i
        if tmpmax >= y[i+1]
            flag = true
            break
        end
    end
    if !flag
        tmpmax = (tmpsum + y[n] - b) / n
    end
    y .= max.(x .- tmpmax, R(0)) # shift, project
    return nothing
end
