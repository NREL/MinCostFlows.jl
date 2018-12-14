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
    imbalances = copy(injections)
    S = Vector{Bool}(undef, n)
    L = Vector{Bool}(undef, n)

    return FlowProblem(
        n, nodesfrom, nodesto, firstfrom, firstto, nextfrom, nextto,
        limits, costs, injections, flows, prices, imbalances, S, L)

end
