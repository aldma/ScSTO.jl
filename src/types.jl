const Maybe{T} = Union{T,Nothing}

"Data for ScSTO problem"
mutable struct ScSTOdata
	# parameters
	x0::Vector{Float64}                         # (augmented) initial state [nx]
    nx::Int                                     # (augmented) state dimension
    N::Int                                      # number of switching intervals
    t0::Float64                                 # initial time
    tf::Float64                                 # desired final time
	dt::Float64									# tf - t0
    swc::Float64                                # switching cost
    Q::Matrix{Float64}                          # state cost matrix [nx x nx]
    E::Matrix{Float64}                          # final state cost matrix [nx x nx]
    U::Matrix{Float64}                       	# inputs sequence [nu x N]
	# functions
    dynam::Function 							# dynamics
    ddynam::Function 							# dynamics derivative
    constr::Maybe{Function} 					# constraints
	dconstr::Maybe{Function} 					# constraints derivative
	pconstr::Maybe{Function} 					# constraints projection
	# discretization
    ngrid::Int                                  # number of fixed grid points
	nvec::Int 									# number of complete grid points
end

"Evaluator for ScSTO problem"
mutable struct ScSTOeval
    # Grid Variables
    tfd::Float64                            	# current final time
    tgrid::Vector{Float64}                      # fixed grid
    tvec::Vector{Float64}                       # complete grid (fixed + sw times)
    tauIdx::Vector{Int}                         # idx of sw times in complete grid
    tgridIdx::Vector{Int}                       # idx of tgrid in complete grid
    dvec::Vector{Float64}              			# complete intervals
    # Shared Data Between Functions
    A::Array{Float64,3}                         # matrices of linearized dynamics
	X::Matrix{Float64}                       	# States at Switching Times
    exm::Array{Float64,3}                    	# Matrix Exponentials
    Phi::Array{Float64,4}                       # Matrix of State Transition Matrices
    M::Array{Float64,3}                         # Integrals over Switching Intervals
    S::Array{Float64,3}                         # S Matrices for each interval
    C::Array{Float64,3}                         # C Matrices for each interval
	jttv::Vector{Float64}
	delta_old_C::Vector{Float64}
	delta_old_S::Vector{Float64}
end

"Reporter for ScSTO problem"
mutable struct ScSTOrepo
	# counters
	ndyna::Int
	nobjf::Int                                  # number of smooth objective evaluations
	ngrad::Int                                  # number of gradient evaluations
	nobjg::Int 									# number of nonsmooth objective evaluations
	nprox::Int 									# number of proximal evaluations
	ncons::Int
	ncjtv::Int
	nproj::Int
	# traces
    objf::Vector{Float64}                      # store smooth objective
    delta::Matrix{Float64}                     # store switching intervals
end

"ScSTO problem"
mutable struct ScSTOptiModel <: AbstractOptiModel
    meta::OptiModelMeta                         # metadata
    data::ScSTOdata 							# problem data
	eval::ScSTOeval                             # evaluator
	repo::ScSTOrepo								# reporter
end

################################################################################
# SHOW, PRINT
################################################################################
import Base.show

show_header(io::IO, prob::ScSTOptiModel) = println(io, typeof(prob))

function show(io::IO, prob::ScSTOptiModel)
	show_header(io, prob)
  	println(io, "Problem name: $(prob.meta.name)")
end
