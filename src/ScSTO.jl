module ScSTO

using OptiMo

export scstoproblem

include("types.jl")
include("interface.jl")
include("l0norm.jl")
include("evaluator.jl")
include("simulation.jl")

end
