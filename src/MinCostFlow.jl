module MinCostFlow

export
    FlowProblem,
    solveflows!,
    updateinjection!,
    updateflowlimit!,
    flows, costs, limits,
    injections, prices

include("FlowProblem.jl")
include("solveflows.jl")
include("update.jl")
include("utils.jl")

end
