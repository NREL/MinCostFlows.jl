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

        #println("Flows: ", flows(fp))
        #println("Prices: ", prices(fp))
        #@assert complementaryslackness(fp)

        # Find a node with positive imbalance
        positiveimbalancenode = firstpositiveimbalance(fp)

        # If no starting node found, problem is either solved or infeasible
        positiveimbalancenode === nothing && break

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
                majoriters_multinode, " multi-node major iterations with ",
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

function decreasereducedcost!(i::Node, ij::Edge{Node}, j::Node, pricechange::Int)

    oldreducedcost = ij.reducedcost
    newreducedcost = oldreducedcost - pricechange # Price change assumed positive on node i
    ij.reducedcost = newreducedcost

    if oldreducedcost === 0 # ij moved from balanced -> active

        # Remove edge from i's balancedfrom adjacency list
        remove!(ij, :prevbalancedfrom, :nextbalancedfrom,
                i, :firstbalancedfrom, :lastbalancedfrom)

        # Remove edge from j's balancedto adjacency list
        remove!(ij, :prevbalancedto, :nextbalancedto,
                j, :firstbalancedto, :lastbalancedto)

        # Add edge to i's activefrom adjacency list
        addend!(ij, :prevactivefrom, :nextactivefrom,
                i, :firstactivefrom, :lastactivefrom)

        # Add edge to j's activeto adjacency list
        addend!(ij, :prevactiveto, :nextactiveto,
                j, :firstactiveto, :lastactiveto)

    elseif newreducedcost === 0 # ij moved from inactive -> balanced

        # Remove edge from i's inactivefrom adjacency list
        remove!(ij, :previnactivefrom, :nextinactivefrom,
                i, :firstinactivefrom, :lastinactivefrom)

        # Remove edge from j's inactiveto adjacency list
        remove!(ij, :previnactiveto, :nextinactiveto,
                j, :firstinactiveto, :lastinactiveto)

        # Add edge to i's balancedfrom adjacency list
        addend!(ij, :prevbalancedfrom, :nextbalancedfrom,
                i, :firstbalancedfrom, :lastbalancedfrom)

        # Add edge to j's balancedto adjacency list
        addend!(ij, :prevbalancedto, :nextbalancedto,
                j, :firstbalancedto, :lastbalancedto)

    end

end

function increasereducedcost!(i::Node, ij::Edge{Node}, j::Node, pricechange::Int)

    oldreducedcost = ij.reducedcost
    newreducedcost = oldreducedcost + pricechange # Price change assumed positive on node j
    ij.reducedcost = newreducedcost

    if oldreducedcost === 0 # ij moved from balanced -> inactive

        # Remove edge from i's balancedfrom adjacency list
        remove!(ij, :prevbalancedfrom, :nextbalancedfrom,
                i, :firstbalancedfrom, :lastbalancedfrom)

        # Remove edge from j's balancedto adjacency list
        remove!(ij, :prevbalancedto, :nextbalancedto,
                j, :firstbalancedto, :lastbalancedto)

        # Add edge to i's inactivefrom adjacency list
        addend!(ij, :previnactivefrom, :nextinactivefrom,
                i, :firstinactivefrom, :lastinactivefrom)

        # Add edge to j's inactiveto adjacency list
        addend!(ij, :previnactiveto, :nextinactiveto,
                j, :firstinactiveto, :lastinactiveto)

    elseif newreducedcost === 0 # ji moved from active -> balanced

        # Remove edge from i's activefrom adjacency list
        remove!(ij, :prevactivefrom, :nextactivefrom,
                i, :firstactivefrom, :lastactivefrom)

        # Remove edge from j's activeto adjacency list
        remove!(ij, :prevactiveto, :nextactiveto,
                j, :firstactiveto, :lastactiveto)

        # Add edge to i's balancedfrom adjacency list
        addend!(ij, :prevbalancedfrom, :nextbalancedfrom,
                i, :firstbalancedfrom, :lastbalancedfrom)

        # Add edge to j's balancedto adjacency list
        addend!(ij, :prevbalancedto, :nextbalancedto,
                j, :firstbalancedto, :lastbalancedto)

    end


end
