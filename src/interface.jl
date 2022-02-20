using LinearAlgebra
using SparseArrays

export warmstart!, initialstate!
export gettau, getdelta

"""
Create ScSTO problem with uniform switching cost
    scstoproblem(x0, dyna, d_dyna, uvec, [constr, d_constr, swc, ngrid, t0, tf, Q, E, tau0ws, name])
"""
function scstoproblem(
    x0::Vector{Float64},                        # initial state [nx]
    dynam::Function,                            # Nonlinear Dynamics
    ddynam::Function,                           # Nonlinear Dynamics Derivative
    U::Matrix{Float64};                         # inputs sequence [nu x N]
    t0::Float64 = 0.0,                          # initial time
    tf::Float64 = 1.0,                          # desired final time
    Q::Maybe{Matrix{Float64}} = nothing,        # state cost matrix [nx x nx]
    E::Maybe{Matrix{Float64}} = nothing,        # final state cost matrix [nx x nx]
    swc::Float64 = 0.0,              # switching cost
    constr::Maybe{Function} = nothing,          # Constraints
    dconstr::Maybe{Function} = nothing,         # Constraints Derivative
    pconstr::Maybe{Function} = nothing,         # Constraints projection
    ncon::Int = 0,                              # number of constraints
    ngrid::Int64 = 100,                         # Number of Linearization points in the fixed grid
    tau0ws::Maybe{Vector{Float64}} = nothing,   # switching times guess [N]
    name::String = "Generic",                   # problem name
)
    @assert ngrid > 1
    @assert tf > t0
    if swc < 0.0
        @warn "negative swc"
    end
    if ncon > 0
        @assert constr !== nothing
        @assert dconstr !== nothing
        @assert pconstr !== nothing
    else
        ncon = 0
    end
    nx = length(x0)                             # state dimension
    N = size(U, 2)                              # number of switching intervals
    if swc === nothing
        swc = 0.0                               # default switching cost
    end
    if Q === nothing
        Q = Diagonal(ones(Float64, nx))         # default cost matrix
    else
        @assert size(Q,1) == nx
        @assert size(Q,2) == nx
    end
    if E === nothing
        E = zeros(Float64, nx, nx)              # default final cost matrix
    else
        @assert size(E,1) == nx
        @assert size(E,2) == nx
    end
    if tau0ws === nothing
        tau0ws = collect(range(t0, stop = tf, length = N + 1))
        tau0ws = tau0ws[2:N]                    # default switching times guess
    else
        @assert length(tau0ws) == N - 1
    end

    # time vectors
    nvec = ngrid + N - 1
    # augmented system state and matrices
    x0 = [x0; 1]
    nx = nx + 1
    spz = Matrix{Float64}(undef, 1, 1)
    spz[1, 1] = 0.0
    spz = sparse(spz)  # Sparse scalar
    Q = Matrix(blockdiag(sparse(Q), spz))
    E = Matrix(blockdiag(sparse(E), spz))

    # problem data
    data = ScSTOdata(
        x0,
        nx,
        N,
        t0,
        tf,
        tf - t0,
        swc,
        Q,
        E,
        U,
        dynam,
        ddynam,
        constr,
        dconstr,
        pconstr,
        ngrid,
        nvec,
    )

    # discretization
    tgrid = collect(range(t0, stop = tf, length = ngrid)) # fixed time grid
    tvec = Vector{Float64}(undef, nvec)         # Complete time grid
    tauIdx = Vector{Int}(undef, N + 1)          # Indeces of switching times in the complete grid
    tgridIdx = Vector{Int}(undef, ngrid)        # Indeces of grid times in the complete grid
    dvec = Vector{Float64}(undef, nvec - 1)     # Intervals of comlete time grid
    # states on complete grid
    X = Matrix{Float64}(undef, nx, nvec)
    X[:, 1] = x0
    # matrices of linearized dynamics
    A = Array{Float64}(undef, nx, nx, nvec - 1)
    exm = Array{Float64}(undef, nx, nx, nvec - 1)
    Phi = Array{Float64}(undef, nx, nx, N + 1, N + 1)
    M = Array{Float64}(undef, nx, nx, nvec - 1)
    S = Array{Float64}(undef, nx, nx, nvec)
    C = Array{Float64}(undef, nx, nx, N)
    delta_old_C = Vector{Float64}(undef, N) # memory, to avoid recomputing
    delta_old_S = Vector{Float64}(undef, N) # memory, to avoid recomputing
    jttv = Vector{Float64}(undef, N-1)
    # evaluator
    eval = ScSTOeval(
        tf,
        tgrid,
        tvec,
        tauIdx,
        tgridIdx,
        dvec,
        A,
        X,
        exm,
        Phi,
        M,
        S,
        C,
        jttv,
        delta_old_C,
        delta_old_S,
    )

    # initial switching intervals
    delta0ws = tau2delta(tau0ws, t0, tf)
    # metadata
    meta = ScSTOmeta(N, ncon, x0=delta0ws, name=name)

    # counters
    ndyna = 0
    nobjf = 0
    ngrad = 0
    nobjg = 0
    nprox = 0
    ncons = 0
    ncjtv = 0
    nproj = 0
    # traces
    objf = Vector{Float64}(undef, 0)
    delta = Matrix{Float64}(undef, N, 0)
    # reporter
    repo = ScSTOrepo(ndyna, nobjf,ngrad, nobjg, nprox, ncons, ncjtv, nproj, objf, delta)

    return ScSTOModel(meta, data, eval, repo)
end

"""
    warmstart!( prob, tau )
"""
function warmstart!(p::ScSTOModel, tau::Vector{Float64})
    @assert length(tau) == p.data.N - 1
    delta = getdelta(p, tau)
    p.meta = ScSTOmeta(
        p.data.N,
        p.meta.ncon,
        x0 = delta,
        y0 = p.meta.y0,
        name = p.meta.name,
    )
    p.repo.objf = Vector{Float64}(undef, 0)
    p.repo.delta = Matrix{Float64}(undef, p.data.N, 0)
end

"""
    initialstate!( prob, x0 )
"""
function initialstate!(p::ScSTOModel, x0::Vector{Float64})
    @assert length(x0) == p.data.nx - 1
    p.data.x0 = [x0; 1]
end

"getdelta(prob, tau)"
function getdelta(p::ScSTOModel, tau::Vector{Float64})
    return tau2delta(tau, p.data.t0, p.data.tf)
end

"gettau(prob, delta)"
function gettau(p::ScSTOModel, delta::Vector{Float64})
    return delta2tau(delta, p.data.t0)
end
