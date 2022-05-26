using ScSTO
using Bazinga

# smooth cost function
struct ScSTOSmoothCost <: ProximableFunction
    model::ScSTOModel
end
function (f::ScSTOSmoothCost)(x)
    return ScSTO.obj(f.model, x)
end
function Bazinga.gradient!(dfx, f::ScSTOSmoothCost, x)
    return ScSTO.objgrad!(f.model, x, dfx)
end

# nonsmooth cost function
struct ScSTONonsmoothCost <: ProximableFunction
    model::ScSTOModel
end
function Bazinga.prox!(z, f::ScSTONonsmoothCost, x, gamma)
    return ScSTO.objprox!(f.model, x, gamma, z)
end

# constraints
struct ScSTOConstraint <: SmoothFunction
    model::ScSTOModel
end
function Bazinga.eval!(cx, f::ScSTOConstraint, x)
    ScSTO.cons!(f.model, x, cx)
    return nothing
end
function Bazinga.jtprod!(jtv, f::ScSTOConstraint, x, v)
    ScSTO.jtprod!(f.model, x, v, jtv)
    return nothing
end

# constraint set
struct ScSTOSet <: ClosedSet
    model::ScSTOModel
end
function Bazinga.proj!(z, f::ScSTOSet, cx)
    ScSTO.proj!(f.model, cx, z)
    return nothing
end
