push!(LOAD_PATH,"/home/alberto/Documents/optimo/src");
push!(LOAD_PATH,"/home/alberto/Documents/bazinga/src");
push!(LOAD_PATH,"/home/alberto/Documents/scsto/src/");

using OptiMo, Bazinga
using ScSTO
using Printf, PyPlot, PyCall

# time interval
t0 = 0.0
tf = 12.0

# control sequence [nu x N]
uvec = Matrix( [repeat([0.0; 1.0], 4, 1); 0.0]' )
#uvec = Matrix( [repeat([0.0; 0.5; 1.0; 0.5], 8, 1); 0.0]' )
nu = size(uvec, 1)
N = size(uvec, 2)

# cost matrix
C = Matrix( [1.0 0.0 -1.0 0.0; 0.0 1.0 0.0  -1.0] )
Q = Matrix( C' * C)

# initial state
x0 = [0.5; 0.7; 1; 1] # nx x 1
nx = length(x0)

# system dynamics
function dynam( x::Vector{Float64}, u::Vector{Float64} )
    n = length(x)
    f = zeros(n)
    f[1] = x[1] - x[1]*x[2] - 0.4*x[1]*u[1]
    f[2] = -x[2] + x[1]*x[2] - 0.2*x[2]*u[1]
    return f
end
function d_dynam( x::Vector{Float64}, u::Vector{Float64} )
    df = [1.0-x[2]-0.4*u[1]       -x[1]                   0    0;
          x[2]                     -1+x[1]-0.2*u[1]       0    0;
          0                        0                      0    0;
          0                        0                      0    0]
    return df
end

# number of points on the fixed grid
ngrid = 100

# switching cost
swc_grid = [0.0; 0.1; 1.0]
#swc_grid = [0.0; 0.0001; 0.0005; 0.005]
ng = length(swc_grid)

# solver
solver=Bazinga.ZEROFPR( tol_optim=1e-4,
                        max_iter=100,
                        verbose=true,
                        freq=10 )

# allocation
nt = 1000
objective = Array{Float64}(undef, ng)
cputime = Array{Float64}(undef, ng)
swdelta = Array{Float64}(undef, N, ng)
swtau = Array{Float64}(undef, N-1, ng)
swtf = Array{Float64}(undef, ng)
xsim = Array{Float64}(undef, nx, nt, ng)
usim = Array{Float64}(undef, nu, nt, ng)

for k in 1:ng

    swc = swc_grid[k]

    # STO problem with switching costs
    prob = scstoproblem(x0, dynam, d_dynam, uvec, ngrid=ngrid, t0=t0, tf=tf, Q=Q, swc=swc);

    # solve!
    out = solver( prob )
    print( out )

    cputime[k] = out.time
    swdelta[:,k] .= out.x
    swtau[:,k], swtf[k] = gettau(prob, swdelta[:,k])

    xsim[:,:,k], _, objective[k], t = simulate(prob, swtau[:,k], length=nt)
    usim[:,:,k], _ = simulateinput(prob, swtau[:,k], t)

end

t = collect( range(t0, stop=tf, length=nt) )

figure()
subplot(3,1,1)
for k in 1:ng
    plot(t, xsim[1,:,k])
end
ylim(0, 1.75)
yticks([0; 1; ])
ylabel(L"x_1")

subplot(3,1,2)
for k in 1:ng
    plot(t, xsim[2,:,k])
end
ylabel(L"x_2")
ylim(0, 1.75)
yticks([0; 1; ])

subplot(3,1,3)
for k in 1:ng
    plot(t, usim[1,:,k])
end
ylim(-0.2, 1.2)
yticks([0; 1])
ylabel(L"u")
xlabel(L"$\mathrm{Time}\; [s]$")

gcf()

swnnz = similar(swc_grid)
for k in 1:ng
    swnnz[k] = sum(swdelta[:,k] .> 0)
end
swobjective = swc_grid .* swnnz
swobjective .+= objective


#=
swc2 = swc_grid[ng]
prob2 = scstoproblem(x0, dynam, d_dynam, uvec, ngrid=ngrid, t0=t0, tf=tf, Q=Q, swc=swc2);
out2 = solver( prob2 )
tau2, _ = gettau(prob2, out2.x)
xsim2, _, _, t2 = simulate(prob2, tau2, length=nt)
usim2, _ = simulateinput(prob2, tau2, t2)

figure()
subplot(3,1,1)
plot(t2, xsim2[1,:])
ylim(0, 1.75)
yticks([0; 1; ])
ylabel(L"x_1")
subplot(3,1,2)
plot(t2, xsim2[2,:])
ylabel(L"x_2")
ylim(0, 1.75)
yticks([0; 1; ])
subplot(3,1,3)
plot(t2, usim2[1,:])
ylim(-0.2, 1.2)
yticks([0; 1])
ylabel(L"u")
xlabel(L"$\mathrm{Time}\; [s]$")

figure()
plot(prob2.repo.objf)
ylabel(L"objective f")
xlabel(L"call")
=#



return nothing
