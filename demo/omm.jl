push!(LOAD_PATH,"/home/albertodm/Documents/OptiMo.jl/src");
push!(LOAD_PATH,"/home/albertodm/Documents/Bazinga.jl/src");
push!(LOAD_PATH,"/home/albertodm/Documents/ScSTO.jl/src/");

using OptiMo, Bazinga
using ScSTO
using Printf, PyPlot, PyCall

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
uvec = Matrix( [repeat([0.0; 1.0; 0.0; 1.0; 0.0; 2.0; 0.0; 1.0], 2, 1); 0.0]' )
nu = size(uvec, 1)
N = size(uvec, 2)

# cost matrices
Q = Matrix( [0.0 0.0; 0.0 0.0] )
E = Matrix( [0.0 0.0; 0.0 -1.0] )

# system dynamics
function dynam( x::Vector{Float64}, u::Vector{Float64} )
    #k = [1.0; 5.0; 1.0; 1.0; 5.0; 4.0; 10.0; 5.0]
    k = [1.0; 3.0; 1.0; 4.0; 2.5; 1.5; 10.0; 2.0]
    if u[1] < 0.5 # no maintenance
        f = [ -k[1]*x[1]; k[2]*x[1] - k[3] ]
    elseif u[1] < 1.5 # minor maintenance
        f = [ k[4]*(1 - x[1]); k[5]*x[1] - k[6] ]
    else # major maintenance
        f = [ k[7]*(1 - x[1]^2); -k[8] ]
    end
    return f
end
function d_dynam( x::Vector{Float64}, u::Vector{Float64} )
    k = [1.0; 3.0; 1.0; 4.0; 2.5; 1.5; 10.0; 2.0]
    if u[1] < 0.5 # no maintenance
        df = [ -k[1]  0 ;
                k[2]  0 ]
    elseif u[1] < 1.5 # minor maintenance
        df = [ -k[4]  0 ;
                k[5]  0 ]
    else # major maintenance
        df = [ -2*k[7]*x[1]  0 ;
                0            0 ]
    end
    return df
end

# switching cost
swc = 0.0 # 0.1

################################################################################
# fixed grid points
ngrid = 200

# STO problem with switching costs
prob = scstoproblem(x0, dynam, d_dynam, uvec, ngrid=ngrid, t0=t0, tf=tf, Q=Q, E=E, swc=swc);

# solve
solver = Bazinga.ZEROFPR( tol_optim=1e-4,
                          max_iter=100,
                          verbose=true,
                          freq=10,
                          gamma_min=1e-16 )
out = solver( prob )
print( out )

################################################################################
# switching intervals
delta0 = copy(prob.meta.x0)
swdelta = out.x
# switching times
tau0, _ = gettau(prob, delta0)
swtau, swtf = gettau(prob, swdelta)
# simulation
nt = 1000
xsim0, _, objective0, tsim = simulate(prob, tau0, length=nt)
usim0, _ = simulateinput(prob, tau0, tsim)
xsim, _, objective, _ = simulate(prob, swtau, tsim)
usim, _ = simulateinput(prob, swtau, tsim)

################################################################################
# plotting
figure()
subplot(3,1,1)
plot(tsim, xsim0[1,:],linestyle="--")#,linewidth=1)
plot(tsim, xsim[1,:])
ylim(0, 1.1)
yticks([0; 1; ])
ylabel(L"x")
subplot(3,1,2)
plot(tsim, xsim0[2,:],linestyle="--")#,linewidth=1)
plot(tsim, xsim[2,:])
ylabel(L"y")
subplot(3,1,3)
plot(tsim, usim0[1,:],linestyle="--")#,linewidth=1)
plot(tsim, usim[1,:])
ylim(-0.2, 2.2)
yticks([0; 1; 2])
ylabel(L"u")
xlabel(L"t")

gcf()

#=
figure()
plot( prob.repo.objf )
ylabel(L"objf")
xlabel(L"call")
=#


return nothing
