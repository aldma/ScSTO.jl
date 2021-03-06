using OptiMo, Bazinga, ScSTO
using Printf, PyPlot, PyCall
using DelimitedFiles

foldername = "/home/albertodm/Documents/ScSTO.jl/demo/"

"""
Maintenance optimization problem

    [SAL14]     Sun, Aw, Loxton, Teo, "An optimal machine maintenance problem with
    probabilistic state constraints" (2014)

    state : machine state `x` and cumulated profit `y`

 3-dynamics
 (i) production : full output, no maintenance
 (ii) minor maintenance : reduced output, slow, cheap maintenance, minor improvement
 (iii) major maintenance : no output, fast maintenance, major improvement
"""
# time interval
t0 = 0.0
tf = 1.0

# initial state
x0 = [1.0; 0.0]
nx = length(x0)

# control sequence
uvec = Matrix(repeat([0.0; 1.0; 0.0; 1.0; 0.0; 2.0], 3, 1)')
uvec = uvec[:, 1:end-1]
nu = size(uvec, 1)
N = size(uvec, 2)

# cost matrices
Q = Matrix([0.0 0.0; 0.0 0.0])
E = Matrix([0.0 0.0; 0.0 -1.0])

# system dynamics
function dynam(x::Vector{Float64}, u::Vector{Float64})
    k = [2.0; 25.0; 1.0; 1.0; 8.0; 2.0; 50.0; 40.0]
    if u[1] < 0.5 # no maintenance
        f = [-k[1] * x[1]; k[2] * x[1]^2 - k[3]]
    elseif u[1] < 1.5 # minor maintenance
        f = [k[4] * (1 - x[1]); k[5] * x[1]^2 - k[6]]
    else # major maintenance
        f = [k[7] * (1 - x[1]^2); -k[8]]
    end
    return f
end
function d_dynam(x::Vector{Float64}, u::Vector{Float64})
    k = [2.0; 25.0; 1.0; 1.0; 8.0; 2.0; 50.0; 40.0]
    if u[1] < 0.5 # no maintenance
        df = [
            -k[1] 0
            2 * k[2] * x[1] 0
        ]
    elseif u[1] < 1.5 # minor maintenance
        df = [
            -k[4] 0
            2 * k[5] * x[1] 0
        ]
    else # major maintenance
        df = [
            -2 * k[7] * x[1] 0
            0 0
        ]
    end
    return df
end

# fixed grid points
ngrid = 400

# STO problem with switching costs
prob = scstoproblem(
    x0,
    dynam,
    d_dynam,
    uvec,
    ngrid = ngrid,
    t0 = t0,
    tf = tf,
    Q = Q,
    E = E,
    swc = 0.0,
);

# solve
solver = Bazinga.ZEROFPR(
    max_iter = 100,
    verbose = true,
    gamma_min = 1e-16,
)
out = solver(prob)
print(out)

# switching intervals
delta0 = copy(prob.meta.x0)
swdelta = out.x
# switching times
tau0, _ = gettau(prob, delta0)
swtau, swtf = gettau(prob, swdelta)
# simulation
nt = 1000
xsim0, _, objective0, tsim = simulate(prob, tau0, length = nt)
usim0, _ = simulateinput(prob, tau0, tsim)
xsim, _, objective, _ = simulate(prob, swtau, tsim)
usim, _ = simulateinput(prob, swtau, tsim)
objf_0 = prob.repo.objf

tmpmat = Array{Float64}(undef, nt, 1+2*3)
tmpmat[:,1] = tsim
tmpmat[:,2] = xsim0[1, :]
tmpmat[:,3] = xsim0[2, :]
tmpmat[:,4] = usim0[1, :]
tmpmat[:,5] = xsim[1, :]
tmpmat[:,6] = xsim[2, :]
tmpmat[:,7] = usim[1, :]

writedlm( foldername * "omm_init_traj.csv",  tmpmat, ',')

#=
figure()
subplot(3, 1, 1)
plot(tsim, xsim0[1, :], linestyle = "--")
plot(tsim, xsim[1, :])
ylim(0, 1.1)
xlim(0, 1)
yticks([0; 1])
ylabel(L"s")
subplot(3, 1, 2)
plot(tsim, xsim0[2, :], linestyle = "--")
plot(tsim, xsim[2, :])
ylabel(L"p")
xlim(0, 1)
subplot(3, 1, 3)
plot(tsim, usim0[1, :], linestyle = "--")
plot(tsim, usim[1, :])
ylim(-0.2, 2.2)
xlim(0, 1)
yticks([0; 1; 2])
ylabel(L"u")
xlabel(L"t")
gcf()
savefig(foldername * "ScSTO.jl/demo/data/omm_traj.pdf")
=#

################################################################################
# with switching cost
################################################################################

swc_grid = [25.0; 30.0] # 20, 25, 30
ng = length(swc_grid)

objective_swc = Array{Float64}(undef, ng)
swdelta_swc = Array{Float64}(undef, N, ng)
swtau_swc = Array{Float64}(undef, N - 1, ng)
xsim_swc = Array{Float64}(undef, nx, nt, ng)
usim_swc = Array{Float64}(undef, nu, nt, ng)
objf_swc = Array{Any}(undef, ng)

for k = 1:ng
    swc = swc_grid[k]

    # STO problem with switching costs
    prob_swc = scstoproblem(
        x0,
        dynam,
        d_dynam,
        uvec,
        ngrid = ngrid,
        t0 = t0,
        tf = tf,
        Q = Q,
        E = E,
        swc = swc,
        tau0ws = swtau,
    )

    out_swc = solver(prob_swc)
    print(out_swc)

    swdelta_swc[:, k] = out_swc.x
    swtau_swc[:, k], _ = gettau(prob_swc, swdelta_swc[:, k])
    xsim_swc[:, :, k], _, objective_swc[k], _ = simulate(prob_swc, swtau_swc[:, k], tsim)
    usim_swc[:, :, k], _ = simulateinput(prob_swc, swtau_swc[:, k], tsim)
    objf_swc[k] = prob_swc.repo.objf
end

tmpmat = Array{Float64}(undef, nt, 1+3*(ng+1))
tmpmat[:,1] = t
tmpmat[:,2] = xsim[1, :]
tmpmat[:,3] = xsim[2, :]
tmpmat[:,4] = usim[1, :]
for k = 1:ng
    tmpmat[:,3*k+2] = xsim_swc[1, :, k]
    tmpmat[:,3*k+3] = xsim_swc[2, :, k]
    tmpmat[:,3*k+4] = usim_swc[1, :, k]
end
writedlm( foldername * "omm_swc_traj.csv",  tmpmat, ',')

#=
figure()
subplot(3, 1, 1)
plot(tsim, xsim[1, :], linestyle = "--")
for k = 1:ng
    plot(tsim, xsim_swc[1, :, k])
end
ylim(0, 1.1)
xlim(0, 1)
yticks([0; 1])
ylabel(L"s")
subplot(3, 1, 2)
plot(tsim, xsim[2, :], linestyle = "--")
for k = 1:ng
    plot(tsim, xsim_swc[2, :, k])
end
ylabel(L"p")
xlim(0, 1)
subplot(3, 1, 3)
plot(tsim, usim[1, :], linestyle = "--", label = "σ = 0.0")
for k = 1:ng
    plot(tsim, usim_swc[1, :, k], label = "σ = $(swc_grid[k])")
end
#legend(loc = "upper right")
legend(loc = "upper center", bbox_to_anchor=(0.5, -0.2), ncol=5)
ylim(-0.2, 2.2)
xlim(0, 1)
yticks([0; 1; 2])
ylabel(L"u")
xlabel(L"t")
gcf()
savefig(foldername * "ScSTO.jl/demo/data/omm_swc_traj.pdf")
=#

#=
figure()
semilogy(objf_0 .- minimum(objf_0))
for k in 1:ng
    semilogy(objf_swc[k] .- minimum(objf_swc[k]))
end
gcf()
savefig(foldername * "ScSTO.jl/demo/omm_objf.pdf")
=#
