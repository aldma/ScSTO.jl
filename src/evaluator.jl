using LinearAlgebra

################################################################################
# Evaluation functions for (unconstrained) ScSTOModel
# obj, grad!, prox!, objprox!
################################################################################
"Evaluate smooth objective"
function obj(p::ScSTOModel, x::Vector{Float64})
    chkfinite(x)
    precompMatJ!(p, x)
    fx = (p.data.x0' * p.eval.S[:, :, 1] * p.data.x0)[1]
    p.repo.nobjf += 1
    p.repo.objf = [p.repo.objf; fx]
    p.repo.delta = [p.repo.delta x]
    chkfinite(fx)
    return fx
end

"Evaluate gradient"
function grad!(p::ScSTOModel, x::Vector{Float64}, dfx::Vector{Float64})
    chkfinite(x)
    precompMatGradJ!(p, x)
    for i = 1:p.data.N
        dfx[i] = (p.eval.X[:, p.eval.tauIdx[i+1]]'*p.eval.C[:, :, i]*p.eval.X[:, p.eval.tauIdx[i+1]])[1]
    end
    p.repo.ngrad += 1
    chkfinite(dfx)
    return nothing
end

"Evaluate smooth objective and gradient"
function objgrad!(p::ScSTOModel, x::Vector{Float64}, dfx::Vector{Float64})
    chkfinite(x)
    precompMatGradJ!(p, x)
    fx = (p.data.x0'*p.eval.S[:, :, 1]*p.data.x0)[1]
    for i = 1:p.data.N
        dfx[i] = (p.eval.X[:, p.eval.tauIdx[i+1]]'*p.eval.C[:, :, i]*p.eval.X[:, p.eval.tauIdx[i+1] ])[1]
    end
    p.repo.nobjf += 1
    p.repo.ngrad += 1
    p.repo.objf = [p.repo.objf; fx]
    p.repo.delta = [p.repo.delta x]
    chkfinite(fx)
    chkfinite(dfx)
    return fx
end

"Constraints"
function cons!(p::ScSTOModel, x::Vector{Float64}, cx::Vector{Float64})
    chkfinite(x)
    tau, _ = gettau(p, x)
    p.data.constr(tau, cx)
    p.repo.ncons += 1
    chkfinite(cx)
    return nothing
end

"Transposed-Jacobian-vector multiplication"
function jtprod!(p::ScSTOModel, x::Vector{Float64}, v::Vector{Float64}, jtv::Vector{Float64})
    chkfinite(x)
    chkfinite(v)
    tau, _ = gettau(p, x)
    # x, jtv [N]
    # tau, jttv [N-1]
    # v [ncon]
    p.data.dconstr(tau, v, p.eval.jttv)
    jtv[p.data.N] = 0.0
    for i in p.data.N-1:-1:1
        jtv[i] = jtv[i+1] + p.eval.jttv[i]
    end
    p.repo.ncjtv += 1
    chkfinite(jtv)
    return nothing
end

"Project constraints"
function proj!(p::ScSTOModel, cx::Vector{Float64}, px::Vector{Float64})
    chkfinite(cx)
    p.data.pconstr(cx, px)
    p.repo.nproj += 1
    chkfinite(px)
    return nothing
end

"Prox operator"
function prox!(p::ScSTOModel, x::Vector{Float64}, a::Float64, z::Vector{Float64})
    chkfinite(x)
    chkfinite(a)
    proxl0simplex!(x, a * p.data.swc, p.data.dt, z)
    p.repo.nprox += 1
    chkfinite(z)
    return nothing
end

"Prox operator and objective"
function objprox!(p::ScSTOModel, x::Vector{Float64}, a::Float64, z::Vector{Float64})
    chkfinite(x)
    chkfinite(a)
    gz = proxl0simplex!(x, a * p.data.swc, p.data.dt, z)
    gz *= p.data.swc
    p.repo.nobjg += 1
    p.repo.nprox += 1
    chkfinite(z)
    chkfinite(gz)
    return gz
end

################################################################################
# Functions to precompute matrices
################################################################################
"Precompute matrices for cost function"
function precompMatJ!(p::ScSTOModel, x::Vector{Float64})
    chkfinite(x)
    if p.eval.delta_old_S != x
        copyto!(p.eval.delta_old_S, x)
        propagateDynamics!(p, x)
        computeS!(p.eval, p.data.E, p.data.nvec)
    end
    return nothing
end

"Precompute matrices for gradient"
function precompMatGradJ!(p::ScSTOModel, x::Vector{Float64})
    chkfinite(x)
    if p.eval.delta_old_C != x
        copyto!(p.eval.delta_old_C, x)
        precompMatJ!(p, x)
        computeC!(p.eval, p.data.Q, p.data.N)
    end
    return nothing
end

"Propagate dynamics"
function propagateDynamics!(p::ScSTOModel, x::Vector{Float64})
    chkfinite(x)
    p.repo.ndyna += 1
    # Get positive delta
    x = max.(x, 0.0)
    # Get switching times from delta
    tau, p.eval.tfd = delta2tau(x, p.data.t0)
    # Create grid from t0 to tfd
    #p.eval.tgrid[end] = p.eval.tfd
    #p.eval.tgrid = min.(p.eval.tgrid, p.eval.tgrid[end])
    p.eval.tgrid .= collect(range(p.data.t0, stop = p.eval.tfd, length=p.data.ngrid))
    chkfinite(p.eval.tgrid)
    # Fit switching times withing the grid
    p.eval.tvec, p.eval.tauIdx, p.eval.tgridIdx = mergeSortFindIndex(p.eval.tgrid, tau)
    p.eval.dvec = diff(p.eval.tvec)
    chkfinite(p.eval.dvec)
    # index for current input
    uIdx = 1
    # Compute Matrix Exponentials
    n = p.data.nx
    for i = 1:p.data.nvec-1
        # Verify which U input applies
        if uIdx < p.data.N && i >= p.eval.tauIdx[uIdx+1]
            uIdx += 1
        end
        # Linearize Dynamics
        chkfinite(p.eval.X[1:end-1, i])
        p.eval.A[:, :, i] = linearizeDyn(
            p.data.dynam,
            p.data.ddynam,
            p.eval.X[1:end-1, i],
            p.data.U[:, uIdx],
        )
        chkfinite(p.eval.A[:, :, i])
        # Compute Matrix Exponential
        tm = Matrix( [ - p.eval.A[:, :, i]' p.data.Q
                       zeros(n, n) p.eval.A[:, :, i] ] * p.eval.dvec[i] )
        chkfinite(tm)
        tm .= exp(tm)
        chkfinite(tm)
        # Assign \mathcal{E}_i
        p.eval.exm[:, :, i] .= tm[n+1:2*n, n+1:2*n]
        # Assign M_k
        p.eval.M[:, :, i] .= tm[n+1:2*n, n+1:2*n]' * tm[1:n, n+1:2*n]
        # Compute state at the next grid point
        p.eval.X[:, i+1] .= p.eval.exm[:, :, i] * p.eval.X[:, i]
    end
    return nothing
end

################################################################################
# Lower level functions
################################################################################
"Low level function for computing S matrices"
function computeS!(e::ScSTOeval, E::Matrix{Float64}, nvec::Int)
    e.S[:, :, nvec] .= E
    for i = nvec-1:-1:1
        e.S[:, :, i] .=
            e.M[:, :, i] .+ e.exm[:, :, i]' * e.S[:, :, i+1] * e.exm[:, :, i]
    end
end

"Low level function for computing C matrices"
function computeC!(e::ScSTOeval, Q::Matrix{Float64}, N::Int)
    for i = 1:N
        e.C[:, :, i] .=
            Q .+
            e.A[:, :, e.tauIdx[i+1]-1]' * e.S[:, :, e.tauIdx[i+1]] .+
            e.S[:, :, e.tauIdx[i+1]] * e.A[:, :, e.tauIdx[i+1]-1]
    end
end

"""
Linearize Dynamics
    A = linearizeDyn(dynam, ddynam, x, u)
"""
function linearizeDyn(dynam::Function, ddynam::Function, x::Vector{Float64}, u::Vector{Float64})
    chkfinite(x)
    chkfinite(u)
    f = dynam(x, u)
    df = ddynam(x, u)
    A = [df f - df * x; zeros(1, length(x) + 1)]
    chkfinite(A)
    return A
end

"""
Merge, sort, and find indices of grid and switching times
    tvec, tauIdx, tgridIdx = mergeSortFindIndex(tgrid, tau)
"""
function mergeSortFindIndex(tgrid::Vector{Float64}, tau::Vector{Float64})
    chkfinite(tgrid)
    chkfinite(tau)
    ngrid = length(tgrid)
    ntau = length(tau)
    tauIdx = Array{Int}(undef, ntau + 2)
    tauIdx[1] = 1
    tauIdx[end] = ntau + ngrid
    tgridIdx = Array{Int}(undef, ngrid)
    tgridIdx[1] = 1
    tgridIdx[end] = ntau + ngrid
    # Create merged and sorted time vector with grid and switching times
    ttmp = vcat(tgrid, tau)     # Concatenate grid and tau vector
    tidxtmp = sortperm(ttmp)    # Find permutation vector to sort ttemp
    tvec = ttmp[tidxtmp]        # Create full sorted tvec
    # Create index of the tau vector elements inside tvec
    for i = 1:ntau
        tauIdx[i+1] = findfirst(x -> x == ngrid + i, tidxtmp)
    end
    # Create index of the tgrid vector elements inside tvec
    for i = 1:ngrid
        tgridIdx[i] = findfirst(x -> x == i, tidxtmp)
    end
    chkfinite(tvec)
    chkfinite(tauIdx)
    chkfinite(tgridIdx)
    return tvec, tauIdx, tgridIdx
end

"""
    delta = tau2delta(tau, t0, tf)
converts from switching times to switching intervals.
"""
function tau2delta(tau::Vector{Float64}, t0::Float64, tf::Float64)
    chkfinite(tau)
    delta = diff([t0; tau; tf])
    chkfinite(delta)
    return delta
end

"""
    tau, tfd = delta2tau(delta, t0)
converts from switching intervals to switching times.
"""
function delta2tau(delta::Vector{Float64}, t0::Float64)
    chkfinite(delta)
    tau = Vector{Float64}(undef, length(delta) - 1)
    tau[1] = t0 + delta[1]
    for i = 2:length(tau)
        tau[i] = tau[i-1] + delta[i]
    end
    tfd = tau[end] + delta[end]
    chkfinite(tau)
    chkfinite(tfd)
    return tau, tfd
end


function chkfinite(x::R) where {R <: Real}
    if !isfinite(x)
        throw(ArgumentError("scalar contains Infs or NaNs"))
    end
    return true
end

function chkfinite(x::Vector{R}) where {R <: Real}
    for e in x
        if !isfinite(e)
            throw(ArgumentError("vector contains Infs or NaNs"))
        end
    end
    return true
end

function chkfinite(x::Matrix{R}) where {R <: Real}
    for e in x
        if !isfinite(e)
            throw(ArgumentError("matrix contains Infs or NaNs"))
        end
    end
    return true
end

function chkfinite(x::Array{R}) where {R <: Real}
    for e in x
        if !isfinite(e)
            throw(ArgumentError("array contains Infs or NaNs"))
        end
    end
    return true
end
