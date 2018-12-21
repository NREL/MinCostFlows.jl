defzero(::Nothing) = 0
defzero(x::Int) = x

function edgelist(fp::FlowProblem)

    edgelist = Vector{Pair{Int,Int}}(undef, length(fp.edges))
    for (e, edge) in enumerate(fp.edges)
        i = nodeidx(fp, edge.nodefrom)
        j = nodeidx(fp, edge.nodeto)
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
        i = nodeidx(fp, node)
        print(i, ", ")
        i in L && error("Loop in L!!")
        push!(L, i)
        node.inS && push!(S, i)
        node = node.nextL
    end

    println("]")
    println("S: ", S)

end

function showSinout(fp::FlowProblem)

    n_edges = length(fp.edges)

    print("intoS: [")
    intoS = 0
    edge = fp.firstintoS
    while edge !== nothing

        if edge.previntoS !== nothing
            i = nodeidx(fp, edge.previntoS.nodefrom)
            j = nodeidx(fp, edge.previntoS.nodeto)
            print("($i => $j) <- ")
        end

        i = nodeidx(fp, edge.nodefrom)
        j = nodeidx(fp, edge.nodeto)
        print("($i => $j)")

        if edge.nextintoS !== nothing
            i = nodeidx(fp, edge.nextintoS.nodefrom)
            j = nodeidx(fp, edge.nextintoS.nodeto)
            print(" -> ($i => $j)")
        end
        print(", ")

        intoS += 1
        intoS > n_edges && error("intoS loop!")
        edge = edge.nextintoS

    end
    println("]")

    print("outofS: [")
    outofS = 0
    edge = fp.firstoutofS
    while edge !== nothing

        if edge.prevoutofS !== nothing
            i = nodeidx(fp, edge.prevoutofS.nodefrom)
            j = nodeidx(fp, edge.prevoutofS.nodeto)
            print("($i => $j) <- ")
        end

        i = nodeidx(fp, edge.nodefrom)
        j = nodeidx(fp, edge.nodeto)
        print("($i => $j)")

        if edge.nextoutofS !== nothing
            i = nodeidx(fp, edge.nextoutofS.nodefrom)
            j = nodeidx(fp, edge.nextoutofS.nodeto)
            print(" -> ($i => $j)")
        end
        print(", ")

        outofS += 1
        outofS > n_edges && error("outofS loop!")
        edge = edge.nextoutofS

    end
    println("]")

    return intoS, outofS

end

nodeidx(fp::FlowProblem, node::Node) = findfirst(n -> n === node, fp.nodes)

function complementarityslackness(fp::FlowProblem)

    satisfied = true

    for edge in fp.edges

        # Check condition for active arcs
        if (edge.reducedcost < 0) && (edge.flow != edge.limit)
            i = findfirst(n -> n === edge.nodefrom, fp.nodes)
            j = findfirst(n -> n === edge.nodeto, fp.nodes)
            println("Active arc violation on edge $i => $j")
            satisfied = false

        # Check condition for balanced arcs
        elseif (edge.reducedcost == 0) && (edge.flow < 0 || edge.flow > edge.limit)
            i = findfirst(n -> n === edge.nodefrom, fp.nodes)
            j = findfirst(n -> n === edge.nodeto, fp.nodes)
            println("Balanced arc violation on edge $i => $j")
            satisfied = false

        # Check condition for inactive arcs
        elseif (edge.reducedcost > 0) && (edge.flow != 0)
            i = findfirst(n -> n === edge.nodefrom, fp.nodes)
            j = findfirst(n -> n === edge.nodeto, fp.nodes)
            println("Inactive arc violation on edge $i => $j")
            satisfied = false

        end

    end

    return satisfied

end
