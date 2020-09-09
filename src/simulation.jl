using DiffEqBase, OrdinaryDiffEq

export simulate, simulateinput

################################################################################
# State trajectory
################################################################################
"simulate(data, tau, t)"
function simulate(data::ScSTOdata, tau::Vector{Float64}, t::Vector{Float64})
    Q = data.Q[1:end-1, 1:end-1]
    E = data.E[1:end-1, 1:end-1]
    x0 = data.x0[1:end-1]
    xsim, xswt, cost = simulateCore(data.dynam, tau, x0, Q, E, data.U, t)
    return xsim, xswt, cost, t
end

"simulate(prob, tau, t)"
function simulate(p::ScSTOptiModel, tau::Vector{Float64}, t::Vector{Float64})
    return simulate(p.data, tau, t)
end

"simulate(prob, tau, [length])"
function simulate(p::ScSTOptiModel, tau::Vector{Float64}; length::Int = 1000)
    @assert 1 < length
    t = collect(range(p.data.t0, stop = p.data.tf, length = length))
    return simulate(p, tau, t)
end

################################################################################
# Control trajectory
################################################################################
"simulateinput(data, tau, t)"
function simulateinput(data::ScSTOdata, tau::Vector{Float64}, t::Vector{Float64})
    u = computeinput(tau, data.U, t)
    return u, t
end

"simulateinput(prob, tau, t)"
function simulateinput(p::ScSTOptiModel, tau::Vector{Float64}, t::Vector{Float64})
    return simulateinput(p.data, tau, t)
end

"simulateinput(prob, tau, [length])"
function simulateinput(p::ScSTOptiModel, tau::Vector{Float64}; length::Int = 1000)
    @assert 1 < length
    t = collect(range(p.data.t0, stop = p.data.tf, length = length))
    return simulateinput(p, tau, t)
end

################################################################################
# LOWER LEVEL FUNCTIONS
################################################################################
"""
Core function for simulating evolution of switched nonlinear dynamics
    xsim, xswt, cost = simulateCore(dynam, tau, x0, Q, E, U, t)
"""
function simulateCore(
    dynam::Function,
    tau::Vector{Float64},
    x0::Vector{Float64},
    Q::Matrix{Float64},
    E::Matrix{Float64},
    U::Matrix{Float64},
    t::Vector{Float64},
)
    nx = length(x0)  # Number of States
    N = length(tau)  # Number of switches
    nt = length(t)
    tau = [t[1]; tau; t[end]]  # Extend tau vector to simplify numbering
    xsim = zeros(nx, length(t)) # allocation
    xsim[:, 1] = x0 # initial state
    # Define indeces to determine current switching mode
    tmpInd1 = 1
    tmpInd2 = 1
    xprevSwitch = x0
    # Create Vector of Points at the switching instants
    xswt = zeros(nx, N + 1)
    xswt[:, 1] = x0
    # ODE solver
    odes = OrdinaryDiffEq.Tsit5()
    for i = 1:N+1 # Integrate over all the intervals
        # redefine Dynamic function
        nldyn(x, p, t) = dynam(x, U[:, i])
        while t[tmpInd2] < tau[i+1] && tmpInd2 < nt
            tmpInd2 = tmpInd2 + 1  # Increase time index
        end
        if tmpInd2 > tmpInd1
            # ODE problem
            odep = DiffEqBase.ODEProblem(nldyn, xprevSwitch, (t[tmpInd1], t[tmpInd2]))
            # ODE solution
            if tmpInd2 == tmpInd1 + 1
                xsol = DiffEqBase.solve(odep, odes; save_everystep = false)
                xtmp = [xsol[1] xsol[end]]
            else
                xsol = DiffEqBase.solve(odep, odes; saveat = t[tmpInd1:tmpInd2])
                xtmp = hcat((xsol[i] for i in eachindex(xsol))...)
            end
            xsim[:, tmpInd1:tmpInd2] = xtmp
            # update
            xprevSwitch = xsim[:, tmpInd2]
            tmpInd1 = tmpInd2
            if i < N + 1
                xswt[:, i+1] = xsim[:, tmpInd2]
            end
        end
    end
    # smooth objective
    cost = trapz(t, diag(xsim' * Q * xsim)) + (xsim[:, end]'*E*xsim[:, end])[1]
    return xsim, xswt, cost
end

"""
Compute Actual Inputs from Artificial Ones
    usim = computeinput(tau, U, t)
"""
function computeinput(tau::Vector{Float64}, U::Matrix{Float64}, t::Vector{Float64})
    N = length(tau) + 1
    nu = size(U, 1) # U [nu x (N+1)]
    nt = length(t)
    tauexp = [t[1]; tau; t[end]] # expand [N+2]
    usim = zeros(nu, length(t)) # allocation [nu x nt]
    # Define indeces to span time vector
    tmpInd1 = 1
    tmpInd2 = 1
    for i = 1:N  # Iterate over all intervals
        while t[tmpInd2] < tauexp[i+1] && tmpInd2 < nt
            tmpInd2 += 1
        end
        if tmpInd2 > tmpInd1
            usim[:, tmpInd1:tmpInd2] = repeat(U[:, i], 1, tmpInd2 - tmpInd1 + 1)
            tmpInd1 = tmpInd2
        end
    end
    return usim
end

"""
Trapezoidal integration rule
    trapz(x, y)
"""
function trapz(x::Vector{Float64}, y::Vector{Float64})
    n = length(x)
    @assert length(y) == n
    r = 0.0
    if n == 1
        return r
    else
        for i = 2:n
            r += (x[i] - x[i-1]) * (y[i] + y[i-1])
        end
        r *= 0.5
    end
    return r
end
