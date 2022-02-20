using ScSTO
using Random
using Printf

R = Float64
ntests = 100
nx = 100

for i in 1:ntests
    a = rand(R)
    b = 2.0 * rand(R)
    x = randn(R, nx)
    y = similar(x)

    nnzy = ScSTO.proxl0nonneg!( x, a, y )
    @assert all(y .>= 0.0)
    @assert sum(y .> 0.0) == nnzy

    nnzy = ScSTO.proxl0simplex!( x, a, b, y )
    @assert isapprox(sum(y), b)
    @assert all(y .>= 0.0)
    @assert sum(y .> 0.0) == nnzy

    ScSTO.projsimplex!(x, b, y)
    @assert isapprox(sum(y), b)
    @assert all(y .>= 0.0)

end

@printf "Passed"
