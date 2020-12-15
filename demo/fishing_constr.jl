using OptiMo, Bazinga, ScSTO
using Printf, PyPlot, PyCall

# time interval
t0 = 0.0
tf = 12.0

# control sequence [nu x N]
uvec = Matrix([repeat([0.0; 1.0], 4, 1); 0.0]')
nu = size(uvec, 1)
N = size(uvec, 2)

# cost matrix
C = Matrix([1.0 0.0 -1.0 0.0; 0.0 1.0 0.0 -1.0])
Q = Matrix(C' * C)

# initial state
x0 = [0.5; 0.7; 1; 1]
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

# constraints
# N-1  ->  ncon
# N-1, ncon  ->  N-1
# ncon  ->  ncon
ncon = 0
function constr(tau::Vector{Float64}, c::Vector{Float64})
    c[1] = tau[2]
    return nothing
end
function d_constr(tau::Vector{Float64}, v::Vector{Float64}, jtv::Vector{Float64})
    jtv .= 0.0
    jtv[2] = v[1]
    return nothing
end
function p_constr(c::Vector{Float64}, p::Vector{Float64})
    p[1] = min(c[1], 4.0)
    return nothing
end

# number of points on the fixed grid
ngrid = 100
nt = 1000

# switching cost
swc = 0.0 # 0.01

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
    swc = swc,
    ncon = ncon,
    constr = constr,
    dconstr = d_constr,
    pconstr = p_constr,
);

# solve unconstrained problem
solver = Bazinga.ZEROFPR(tol_optim = 1e-4, max_iter = 100, verbose = true, freq = 10)
out = solver(prob)
print(out)

swdelta = out.x
swtau, swtf = gettau(prob, swdelta)

xsim, _, objective, t = simulate(prob, swtau, length = nt)
usim, _ = simulateinput(prob, swtau, t)





# solve constrained problem
solver2 = Bazinga.ALPX(
    tol_optim = 1e-6,
    tol_cviol = 1e-6,
    max_iter = 5,
    max_sub_iter = 20,
    verbose = true,
)
warmstart!(prob, swtau)
if ncon > 0
    out2 = solver2(prob)
    print(out2)

    swdelta2 = out2.x
    swtau2, _ = gettau(prob, swdelta2)
    xsim2, _, _, _ = simulate(prob, swtau2, t)
    usim2, _ = simulateinput(prob, swtau2, t)
end


# plotting
figure()
subplot(3, 1, 1)
plot(t, xsim[1, :])
if ncon > 0
    plot(t, xsim2[1, :])
end
ylim(0, 1.75)
yticks([0; 1])
ylabel(L"x_1")
subplot(3, 1, 2)
plot(t, xsim[2, :])
if ncon > 0
    plot(t, xsim2[2, :])
end
ylabel(L"x_2")
ylim(0, 1.75)
yticks([0; 1])
subplot(3, 1, 3)
plot(t, usim[1, :])
if ncon > 0
    plot(t, usim2[1, :])
end
ylim(-0.2, 1.2)
yticks([0; 1])
ylabel(L"u")
xlabel(L"$\mathrm{Time}\; [s]$")
gcf()



return nothing
