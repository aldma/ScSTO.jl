push!(LOAD_PATH,"/home/albertodm/Documents/OptiMo.jl/src");
push!(LOAD_PATH,"/home/albertodm/Documents/Bazinga.jl/src");
push!(LOAD_PATH,"/home/albertodm/Documents/ScSTO.jl/src/");

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
uvec = Matrix( [repeat([0.0; 1.0], 5, 1); 0.0]' )
nu = size(uvec, 1)
N = size(uvec, 2)

# cost matrix
C = Matrix( [1.0 0.0 -1.0 0.0; 0.0 1.0 0.0  -1.0] )
Q = Matrix( C' * C)

# initial state
x0 = [0.5; 0.7; 1; 1] # nx x 1
nx = length(x0)

# system dynamics
# nx, nu  ->  nx
# nx, nu  -> nx x nx
function dynam( x::Vector{Float64}, u::Vector{Float64} )
    n = length(x)
    f = zeros(n)
    f[1] = x[1] - x[1]*x[2] - 0.4*x[1]*u[1]
    f[2] = -x[2] + x[1]*x[2] - 0.2*x[2]*u[1]
    return f
end
function d_dynam( x::Vector{Float64}, u::Vector{Float64} )
    n = length(x)
    df = zeros(n,n)
    df[1,1] = 1.0 - x[2] - 0.4*u[1]
    df[1,2] = -x[1]
    df[2,1] = x[2]
    df[2,2] = -1.0 + x[1] - 0.2*u[1]
    return df
end

# number of points on the fixed grid
ngrid = 100

# switching cost
swc_grid = [0.0; 0.1; 1.0]
col_grid = [:steelblue, :orangered, :green]
ng = length(swc_grid)

# solver
solver=Bazinga.ZEROFPR( tol_optim=1e-4,
                        max_iter=50,
                        verbose=true,
                        freq=10 )

# allocation
nt = 1000
objective = Array{Float64}(undef, ng, 2)
cputime = Array{Float64}(undef, ng, 2)
swdelta = Array{Float64}(undef, N, ng, 2)
swtau = Array{Float64}(undef, N-1, ng, 2)
swtf = Array{Float64}(undef, ng, 2)
xsim = Array{Float64}(undef, nx, nt, ng, 2)
usim = Array{Float64}(undef, nu, nt, ng, 2)
repo_objf = Array{Any}(undef, ng, 2)

for k in 1:ng

    swc = swc_grid[k]

    # STO problem with switching costs
    prob = scstoproblem(x0, dynam, d_dynam, uvec, ngrid=ngrid, t0=t0, tf=tf, Q=Q, swc=swc);

    # solve!
    out = solver( prob )
    print( out )

    cputime[k,1] = out.time
    swdelta[:,k,1] .= out.x
    swtau[:,k,1], swtf[k,1] = gettau(prob, swdelta[:,k,1])

    xsim[:,:,k,1], _, objective[k,1], t = simulate(prob, swtau[:,k,1], length=nt)
    usim[:,:,k,1], _ = simulateinput(prob, swtau[:,k,1], t)

    repo_objf[k,1] = prob.repo.objf

end

t = collect( range(t0, stop=tf, length=nt) )

figure()
subplot(3,1,1)
for k in 1:ng
    plot(t, xsim[1,:,k,1],c=col_grid[k])
end
ylim(0, 1.75)
xlim(0, 12)
yticks([0; 1; ])
ylabel(L"x_1")
subplot(3,1,2)
for k in 1:ng
    plot(t, xsim[2,:,k,1],c=col_grid[k])
end
ylabel(L"x_2")
ylim(0, 1.75)
xlim(0, 12)
yticks([0; 1; ])
subplot(3,1,3)
for k in 1:ng
    plot(t, usim[1,:,k,1],label="$(swc_grid[k])",c=col_grid[k])
end
legend()
ylim(-0.2, 1.2)
yticks([0; 1])
xlim(0, 12)
ylabel(L"u")
xlabel(L"$\mathrm{Time}\; [s]$")
gcf()

#savefig("/home/alberto/Documents/ScSTO.jl/demo/fishing_unc_traj.pdf")

objfmin = minimum(repo_objf[1,1])
for k in 2:ng
    global objfmin = min(objfmin, minimum(repo_objf[k,1]) )
end

figure()
for k in 1:ng
    semilogy(repo_objf[k,1] .- objfmin)
end
gcf()

#figure()
#for k in 1:ng
#    semilogy(repo_objf[k,1] .- minimum(repo_objf[k,1]))
#end
#gcf()

#swnnz = similar(swc_grid)
#for k in 1:ng
#    swnnz[k] = sum(swdelta[:,k,1] .> 0)
#end
#swobjective = swc_grid .* swnnz
#swobjective .+= objective

################################################################################
################################################################################
# constrained STO
################################################################################
################################################################################
# constraints
# N-1  ->  ncon
# N-1, ncon  ->  N-1
# ncon  ->  ncon
ncon = 2
function constr( tau::Vector{Float64}, c::Vector{Float64} )
    c[1] = tau[2]
    c[2] = tau[4]
    return nothing
end
function d_constr( tau::Vector{Float64}, v::Vector{Float64}, jtv::Vector{Float64} )
    jtv .= 0.0
    jtv[2] = v[1]
    jtv[4] = v[2]
    return nothing
end
function p_constr( c::Vector{Float64}, p::Vector{Float64} )
    p[1] = max(1.0, min(c[1], 2.0))
    p[2] = max(3.0, min(c[2], 4.0))
    return nothing
end

solver2 = Bazinga.ALPX( tol_optim=1e-4,
                        tol_cviol=1e-6,
                        max_iter=10,
                        max_sub_iter=10,
                        subsolver=:zerofpr,
                        verbose=true )

for k in 1:ng

    swc = swc_grid[k]

    prob2 = scstoproblem(x0, dynam, d_dynam, uvec, ngrid=ngrid, t0=t0, tf=tf, Q=Q, swc=swc,
                         ncon=ncon, constr=constr, dconstr=d_constr, pconstr=p_constr );

    warmstart!(prob2, swtau[:,k,1])

    out2 = solver2( prob2 )
    print( out2 )

    swdelta[:,k,2] = out2.x
    swtau[:,k,2], _ = gettau(prob2, swdelta[:,k,2])
    xsim[:,:,k,2], _, _, _ = simulate(prob2, swtau[:,k,2], t)
    usim[:,:,k,2], _ = simulateinput(prob2, swtau[:,k,2], t)
end

figure()
subplot(3,1,1)
k = 2
plot(t, xsim[1,:,k,1], ls=:dashed)
plot(t, xsim[1,:,k,2])
ylim(0, 1.75)
xlim(0, 12)
yticks([0; 1; ])
ylabel(L"x_1")
subplot(3,1,2)
plot(t, xsim[2,:,k,1], ls=:dashed)
plot(t, xsim[2,:,k,2])
ylabel(L"x_2")
ylim(0, 1.75)
xlim(0, 12)
yticks([0; 1; ])
subplot(3,1,3)
plot(t, usim[1,:,k,1], ls=:dashed, label="$(swc_grid[k]) unc")
plot(t, usim[1,:,k,2], label="$(swc_grid[k]) con")
legend(loc="upper right",ncol=2)
ylim(-0.2, 1.2)
xlim(0, 12)
yticks([0; 1])
ylabel(L"u")
xlabel(L"$\mathrm{Time}\; [s]$")
gcf()
#savefig("/home/alberto/Documents/ScSTO.jl/demo/fishing_con_traj.pdf")

return nothing
