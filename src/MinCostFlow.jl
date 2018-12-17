module MinCostFlow

export
    FlowProblem,
    solveflows!,
    updateinjection!,
    updateflowlimit!

include("FlowProblem.jl")
include("solveflows.jl")
include("update.jl")
include("utils.jl")

end
