"""
Implements the Relaxation dual-ascent method for solving
min-cost flow problems, as described in Bertsekas (1998)
"""
function solveflows!(fp::FlowProblem; verbose::Bool=false)

    #elist = edgelist(fp)
    starttime = time()

    calculateimbalances!(fp)

    #println(edgelist(fp))
    #println("Costs: ", costs(fp))
    #println("Limits: ", limits(fp))
    #println("Injections: ", injections(fp))

    majoriters = 0
    majoriters_multinode = 0
    minoriters = 0

    while true

        # Find a node with positive imbalance
        positiveimbalancenode = firstpositiveimbalance(fp)

        # If no starting node found, problem is either solved or infeasible
        positiveimbalancenode === nothing && break

        #println("Flows: ", flows(fp))
        #println("Prices: ", prices(fp))
        #@assert complementaryslackness(fp)
        #println("Starting major iteration")

        majoriters += 1

        # First try a single-node relaxation iteration,
        # moving to next iteration if it makes changes
        singlenodeupdate!(fp, positiveimbalancenode) && continue

        # If the single-node iteration doesn't change anything, 
        # run a multi-node relaxation iteration, updating either
        # flows and/or shadow prioces
        majoriters_multinode += 1
        minoriters += multinodeupdate!(fp, positiveimbalancenode)

    end

    endtime = time()
    if verbose
        elapsedtime = endtime - starttime
        majoriters_singlenode = majoriters - majoriters_multinode
        println(majoriters, " major iterations: ",
                majoriters_singlenode, " single-node major iterations, ",
                majoriters_multinode, "multi-node major iterations with ",
                minoriters, " multi-node minor iterations")
        #println("Average ", minoriters/majoriters, " minor iterations per major multi-node iteration")
        println("Solved in ", elapsedtime, "s")
    end

    return fp

end

#TODO: A linked list that gets trimmed down on-the-fly might be faster here

"""
Finds first node in `nodes` with a positive imbalance
"""
function firstpositiveimbalance(fp::FlowProblem)
    for node in fp.nodes
        (node.imbalance > 0) && return node
    end
    return nothing
end

# TODO: Should probably recalculate these on the fly instead (whenever flow is updated)
function calculateimbalances!(fp::FlowProblem)

    for node in fp.nodes
        node.imbalance = node.injection
    end

    for edge in fp.edges
        edge.nodefrom.imbalance -= edge.flow
        edge.nodeto.imbalance += edge.flow
    end

end

# #TODO: Should do this on the fly / with iterative updates instead
function calculateascentgradient!(fp::FlowProblem)

    fp.ascentgradient = 0
    i = fp.firstS
    while i !== nothing

        ij = i.firstbalancedfrom
        while ij !== nothing
            !ij.nodeto.inS && (fp.ascentgradient -= ij.limit)
            ij = ij.nextbalancedfrom
        end

        ij = i.firstactivefrom
        while ij !== nothing
            !ij.nodeto.inS && (fp.ascentgradient -= ij.limit)
            ij = ij.nextactivefrom
        end

        ji = i.firstactiveto
        while ji !== nothing
            !ji.nodefrom.inS && (fp.ascentgradient += ji.limit)
            ji = ji.nextactiveto
        end

        fp.ascentgradient += i.injection

        i = i.nextS

    end

end
