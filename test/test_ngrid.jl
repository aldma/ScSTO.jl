foldername = "/home/alberto/Documents/"
push!(LOAD_PATH, foldername * "OptiMo.jl/src");
push!(LOAD_PATH, foldername * "Bazinga.jl/src");
push!(LOAD_PATH, foldername * "ScSTO.jl/src/");

using OptiMo, Bazinga
using ScSTO
using Printf
using PyPlot, PyCall

################################################################################
################################################################################
# unconstrained STO
################################################################################
################################################################################
# time interval
t0 = 0.0
tf = 12.0

# control sequence [nu x N]
uvec = Matrix([repeat([0.0; 1.0], 5, 1); 0.0]')
nu = size(uvec, 1)
N = size(uvec, 2)

# cost matrix
C = Matrix([1.0 0.0 -1.0 0.0; 0.0 1.0 0.0 -1.0])
Q = Matrix(C' * C)

# initial state
x0 = [0.5; 0.7; 1; 1] # nx x 1
nx = length(x0)

# system dynamics
# nx, nu  ->  nx
# nx, nu  -> nx x nx
function dynam(x::Vector{Float64}, u::Vector{Float64})
    n = length(x)
    f = zeros(n)
    f[1] = x[1] - x[1] * x[2] - 0.4 * x[1] * u[1]
    f[2] = -x[2] + x[1] * x[2] - 0.2 * x[2] * u[1]
    return f
end
function d_dynam(x::Vector{Float64}, u::Vector{Float64})
    n = length(x)
    df = zeros(n, n)
    df[1, 1] = 1.0 - x[2] - 0.4 * u[1]
    df[1, 2] = -x[1]
    df[2, 1] = x[2]
    df[2, 2] = -1.0 + x[1] - 0.2 * u[1]
    return df
end

# switching cost
swc = 0.1

# number of points on the fixed grid
ngrid_grid = [100; 200; 1000]
ng = length(ngrid_grid)

# solver
solver = Bazinga.ZEROFPR( max_iter = 50, verbose = true )

# allocation
nt = 2000;
swdelta = Array{Float64}(undef, N, ng);
swtau = Array{Float64}(undef, N - 1, ng);
xsim = Array{Float64}(undef, nx, nt, ng);
usim = Array{Float64}(undef, nu, nt, ng);

for k = 1:ng

    ngrid = ngrid_grid[k]

    prob = scstoproblem(
        x0,
        dynam,
        d_dynam,
        uvec,
        ngrid = ngrid,
        t0 = t0,
        tf = tf,
        Q = Q,
        swc = swc,
    )

    out = solver(prob)
    print(out)

    swdelta[:, k] .= out.x
    swtau[:, k], _ = gettau(prob, swdelta[:, k])

    xsim[:, :, k], _, _, t = simulate(prob, swtau[:, k], length = nt)
    usim[:, :, k], _ = simulateinput(prob, swtau[:, k], t)

    @printf "\n\n"

end

t = collect(range(t0, stop = tf, length = nt))

figure()
subplot(3, 1, 1)
for k = 1:ng
    plot(t, xsim[1, :, k])
end
ylim(0, 1.75)
xlim(0, 12)
yticks([0; 1])
ylabel(L"x_1")
subplot(3, 1, 2)
for k = 1:ng
    plot(t, xsim[2, :, k])
end
ylabel(L"x_2")
ylim(0, 1.75)
xlim(0, 12)
yticks([0; 1])
subplot(3, 1, 3)
for k = 1:ng
    plot(t, usim[1, :, k], label = "n = $(ngrid_grid[k])")
end
legend()
ylim(-0.2, 1.2)
yticks([0; 1])
xlim(0, 12)
ylabel(L"u")
xlabel(L"$\mathrm{Time}\; [s]$")
gcf()
savefig(foldername * "ScSTO.jl/test/data/fishing_ngrid.pdf")


################################################################################
# make switching times vectors comparable and compare results

using LinearAlgebra

function cleandelta(delta::Vector{Float64})
    # !!! for binary controls only !!!
    # uvec = 0, 1, 0, 1, 0, 1, ...
    deltanew = copy( delta )
    n = length( delta )
    for k in 1:n-2
        if deltanew[k] > 0.0 && deltanew[k+1] == 0.0
            deltanew[k+2] += deltanew[k]
            deltanew[k] = 0.0
        end
    end
    return deltanew
end
function cleandelta(delta::Matrix{Float64})
    deltanew = copy( delta )
    m = size( delta,2 )
    for k in 1:m
        deltanew[:,k] .= cleandelta( delta[:,k] )
    end
    return deltanew
end

swdelta_clean = cleandelta( swdelta )
global swdelta_diff = Array{Float64}(undef, ng-1);
for k in 1:ng-1
    swdelta_diff[k] = norm( swdelta_clean[:,k] - swdelta_clean[:,ng] , Inf )
end
