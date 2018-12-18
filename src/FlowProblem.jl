# TODO: Create a single vector of structs containing edge/node data
#       to avoid getindex'ing all the time

#mutable struct Edge
#    nodefrom::Node
#    nodeto::Node
#    nextfrom::Edge
#    nextto::Edge
#    limit::Int
#    cost::Int
#    reducedcost::Int
#    flow::Int
#    intoS::Bool
#    outofS::Bool
#end

#mutable struct Node
#    firstto::Edge
#    nextto::Edge
#    injection::Int
#    price::Int
#    imbalance::Int
#    inS::Bool
#    inL::Bool
#end

struct FlowProblem
   nodes::Int
   nodesfrom::Vector{Int}
   nodesto::Vector{Int}
   firstfrom::Vector{Int}
   firstto::Vector{Int}
   nextfrom::Vector{Int}
   nextto::Vector{Int}
   limits::Vector{Int}
   costs::Vector{Int}
   injections::Vector{Int}
   flows::Vector{Int}
   shadowprices::Vector{Int}
   imbalances::Vector{Int}
   S::Vector{Bool}
   L::Vector{Bool}
   labels::Vector{Int}
   forwardedges::Vector{Bool}
   backwardedges::Vector{Bool}
   reducedcosts::Vector{Int}
   transS::Vector{Int}
end

function FlowProblem(nodesfrom::Vector{Int}, nodesto::Vector{Int},
                     limits::Vector{Int}, costs::Vector{Int},
                     injections::Vector{Int})

    n = length(injections)
    e = length(nodesfrom)

    @assert length(nodesto) == e
    @assert length(limits) == e
    @assert length(costs) == e
    @assert sum(injections) == 0

    firstfrom = zeros(Int, n)
    firstto = zeros(Int, n)
    nextfrom = Vector{Int}(undef, e)
    nextto = Vector{Int}(undef, e)

    for i in 1:e

        nodefrom = nodesfrom[i]
        nodeto = nodesto[i]

        firstfrom[nodefrom] == 0 && (firstfrom[nodefrom] = i)
        firstto[nodeto] == 0 && (firstto[nodeto] = i)

        nextfrom[i] = defzero(findnext(x -> x == nodefrom, nodesfrom, i+1))
        nextto[i] = defzero(findnext(x -> x == nodeto, nodesto, i+1))

    end

    # Initialize flows and prices to zero:
    # imbalances are therefore just injections
    flows = zeros(Int, e)
    prices = zeros(Int, n)
    imbalances = Vector{Int}(undef, n)

    S = Vector{Bool}(undef, n)
    L = Vector{Bool}(undef, n)
    labels = Vector{Int}(undef, n)
    forwardedges = Vector{Bool}(undef, e)
    backwardedges = Vector{Bool}(undef, e)
    reducedcosts = copy(costs)
    transS = Vector{Int}(undef, e)

    return FlowProblem(
        n, nodesfrom, nodesto, firstfrom, firstto, nextfrom, nextto,
        limits, costs, injections, flows, prices, imbalances,
        S, L, labels, forwardedges, backwardedges, reducedcosts, transS)

end
