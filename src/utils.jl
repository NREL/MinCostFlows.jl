defzero(::Nothing) = 0
defzero(x::Int) = x

function edgelist(fp::FlowProblem)

    edgelist = Vector{Pair{Int,Int}}(undef, length(fp.edges))
    for (e, edge) in enumerate(fp.edges)
        i = findfirst(n -> n === edge.nodefrom, fp.nodes)
        j = findfirst(n -> n === edge.nodeto, fp.nodes)
        edgelist[e] = Pair(i,j)
    end

    return edgelist

end

function showSL(fp::FlowProblem)

    S = Int[]
    L = Int[]
    print("L: [")

    node = fp.firstL
    while node !== nothing
        i = findfirst(n -> n === node, fp.nodes)
        print(i, ", ")
        i in L && error("Loop in L!!")
        push!(L, i)
        node.inS && push!(S, i)
        node = node.nextL
    end

    println("]")
    println("S: ", S)

end

function lengthSinout(fp::FlowProblem)

    n_edges = length(fp.edges)

    intoS = 0
    edge = fp.firstintoS
    while edge !== nothing
        intoS += 1
        intoS > n_edges && error("intoS loop!")
        edge = edge.nextintoS
    end

    outofS = 0
    edge = fp.firstoutofS
    while edge !== nothing
        outofS += 1
        outofS > n_edges && error("outofS loop!")
        edge = edge.nextoutofS
    end

    return intoS, outofS

end

function complementarityslackness(fp::FlowProblem)
    for edge in fp.edges

        # Check condition for active arcs
        if (edge.reducedcost < 0) && (edge.flow != edge.limit)
            println("Active arc violation")
            return false

        # Check condition for balanced arcs
        elseif (edge.reducedcost == 0) && (edge.flow < 0 || edge.flow > edge.limit)
            println("Balanced arc violation")
            return false

        # Check condition for inactive arcs
        elseif (edge.reducedcost > 0) && (edge.flow != 0)
            println("Inactive arc violation")
            return false

        end

    end
    return true
end
