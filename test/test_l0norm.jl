using Random

include("/home/alberto/Documents/scsto/src/l0norm.jl")

R = Float64
ntests = 100
nx = 100

for i in 1:ntests
    a = rand(R)
    b = 2.0 * rand(R)
    x = randn(R, nx)
    y = similar(x)

    nnzy = proxl0nonneg!( x, a, y )
    @assert all(y .>= 0.0)
    @assert sum(y .> 0.0) == nnzy

    nnzy = proxl0simplex!( x, a, b, y )
    @assert isapprox(sum(y), b)
    @assert all(y .>= 0.0)
    @assert sum(y .> 0.0) == nnzy

    projsimplex!(x, b, y)
    @assert isapprox(sum(y), b)
    @assert all(y .>= 0.0)

end

@printf "Passed"
