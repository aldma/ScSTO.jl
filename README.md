# ScSTO.jl
Sparse and Constrained Switching Time Optimization

ScSTO is a modeling tool for switching time optimization (STO) problems in [Julia](https://julialang.org/).
ScSTO builds upon [SwitchTimeOpt.jl](https://github.com/oxfordcontrol/SwitchTimeOpt.jl) and extends it, in that nonlinear switched systems, switching costs, and constraints on the switching times can be easily incorporated.

ScSTO generates an instance of ´´´ScSTOptiModel´´´, which fit into the framework offered by [OptiMo](https://github.com/aldma/OptiMo.jl). Thus, such models can be solved by invoking solvers from [Bazinga](https://github.com/aldma/Bazinga.jl).
