"""
Implements the Relaxation dual-ascent method for solving
min-cost flow problems, as described in Bertsekas (1998)
"""
function solveflows!(fp::FlowProblem; verbose::Bool=true)

    #elist = edgelist(fp)
    starttime = time()

    calculateimbalances!(fp)

    majoriters = 0
    minoriters = 0
    while true

        majoriters += 1

        #println(edgelist(fp))
        #println("Costs: ", costs(fp))
        #println("Limits: ", limits(fp))
        #println("Flows: ", flows(fp))
        #println("Injections: ", injections(fp))
        #println("Prices: ", prices(fp))
        #@assert complementaryslackness(fp)
        #println("Starting major iteration")

        # Find a potential starting node for an augmenting path
        augmentingpathstart = firstpositiveimbalance(fp)

        # No starting node found, problem is either solved or infeasible
        augmentingpathstart === nothing && break

        # Reset the scan sets for this iteration
        resetSL!(fp, augmentingpathstart)
        #showSL(fp)

        # Main iteration: update either flows or shadow prices
        minoriters += update!(fp, augmentingpathstart)

    end

    endtime = time()
    if verbose
        elapsedtime = endtime - starttime
        println(majoriters, " major iterations, ", minoriters, " minor iterations")
        println("Average ", minoriters/majoriters, " minor iterations per major iteration")
        println("Solved in ", elapsedtime, "s")
    end

    return fp

end

function update!(fp::FlowProblem, augmentingpathstart::Node)

    iters = 0
    while true

        iters += 1
        #println("Starting minor iteration")

        # Look for a candidate node i to scan and add to S
        i = augmentS!(fp)
        #showSL(fp)

        # Update prices if it will improve the dual solution
        # or if there are no nodes left to scan
        if i === nothing || dualascendable(fp)
            updateprices!(fp)
            break
        end

        # Label neighbour nodes of i
        augmentingpathend = augmentL!(fp, i)
        #showSL(fp)

        # Didn't find an augmenting path, try adding a different node
        augmentingpathend === nothing && continue

        # Found an augmenting path, augment flows accordingly
        updateflows!(fp, augmentingpathstart, augmentingpathend)
        break

    end

    return iters

end

"""
Finds first node in `nodes` with a positive imbalance
"""
function firstpositiveimbalance(fp::FlowProblem)
    for node in fp.nodes
        (node.imbalance > 0) && return node
    end
    return nothing
end


"""
Empties the set S and reduces the set L to a single element, `j`
"""
function resetSL!(fp::FlowProblem, j::Node)

    #println("Resetting SL...")
    node = fp.firstL
    while node !== nothing
        node.inL = false
        node.inS = false
        node = node.nextL
    end

    # Reset L to {j}
    fp.firstL = j
    j.inL = true
    j.nextL = nothing

    return nothing

end

"""
If S is a subset of L, return an element of L - S after adding that element
to S. If S == L, return `nothing` and make no changes.
"""
function augmentS!(fp::FlowProblem)

    #println("Augmenting S...")
    i = fp.firstL
    while i !== nothing

        if !i.inS

            # i joins S
            i.inS = true

            ij = i.firstbalancedfrom
            while ij !== nothing
                if ij.nodeto.inS
                    fp.ascentgradient += ij.flow
                else
                    fp.ascentgradient -= (ij.limit - ij.flow)
                end
                ij = ij.nextbalancedfrom
            end

            ji = i.firstbalancedto
            while ji !== nothing
                if ji.nodefrom.inS
                    fp.ascentgradient += (ji.limit - ji.flow)
                else
                    fp.ascentgradient -= ji.flow
                end
                ji = ji.nextbalancedfrom
            end

            return i

        end

        i = i.nextL

    end

    return nothing

end

"""
Checks whether the current set S would increase the
dual objective function if prices were updated
"""
function dualascendable(fp::FlowProblem)

    #println("Gradient: ", fp.ascentgradient)
    return fp.ascentgradient > 0

end

"""
Add nodes to L that are adjacent to `i` and satisfy inclusion criteria,
returning one such added node that has a positive imbalance
(or `nothing` if none exist). Subset of Step 2 in Bertsekas.
"""
function augmentL!(fp::FlowProblem, i::N)::Union{N,Nothing} where {N<:Node}

    #println("Augmenting L...")
    negativeimbalancenode = nothing

    # Iterate through balanced edges ji connecting j to i
    ji = i.firstbalancedto
    while ji !== nothing

        j = ji.nodefrom # source node of edge ji

        if (!j.inL && ji.flow > 0) # flow enters L and exceeds lower bound

            # Add j to L
            j.inL = true
            j.nextL = fp.firstL
            fp.firstL = j

            # Save path back to i
            j.augpathprev = ji

            # Check for negative imbalance and save for returning
            j.imbalance < 0 && (negativeimbalancenode = j)

        end

        ji = ji.nextbalancedto # Get next edge

    end

    # Iterate through balanced edges ij connecting i to j
    ij = i.firstbalancedfrom
    while ij !== nothing

        j = ij.nodeto # sink node of edge ij

        if (!j.inL && ij.flow < ij.limit) # flow leaves L and is less than limit

            # Add j to L
            j.inL = true
            j.nextL = fp.firstL
            fp.firstL = j

            # Save path back to i
            j.augpathprev = ij

            # Check for negative imbalance and save for returning
            j.imbalance < 0 && (negativeimbalancenode = j)

        end

        ij = ij.nextbalancedfrom # Get next edge

    end

    return negativeimbalancenode

end

"""
Updates the flows on the augmenting path from `startnode` to `endnode`
(Step 3 in Bertsekas)
"""
function updateflows!(fp::FlowProblem, startnode::Node, endnode::Node)

    #println("Updating flows...")
    delta = min(startnode.imbalance, -endnode.imbalance)

    # First pass, determine value of delta
    # Traverse the augmenting path backwards to startnode
    currentnode = endnode
    while currentnode !== startnode

        # Determine edge direction and move to previous node
        ab = currentnode.augpathprev # Get edge
        a = ab.nodefrom # Get edge to/from nodes
        b = ab.nodeto

        if b === currentnode # ab is a forward edge in the path, move to node a
            delta = min(delta, ab.limit - ab.flow)
            currentnode = a
        else # ab is a backwards edge in the path, move to node b
            delta = min(delta, ab.flow)
            currentnode = b
        end

    end

    # Second pass, adjust flows by delta
    currentnode = endnode
    while currentnode !== startnode

        # Determine edge direction and move to previous node
        ab = currentnode.augpathprev # Get edge
        a = ab.nodefrom # Get edge to/from nodes
        b = ab.nodeto

        if b === currentnode # ab is a forward edge, move to node a
            ab.flow += delta
            currentnode = a
        else # ab is a backward edge, move to node b
            ab.flow -= delta
            currentnode = b
        end

    end

    calculateimbalances!(fp) #TODO: Do this on the fly instead

end

"""
Updates the shadow prices (and potentially flows) of the elements of S
(Step 4 in Bertsekas)
"""
function updateprices!(fp::FlowProblem)

    function printedge(edge::Edge)
        a = nodeidx(fp, edge.nodefrom)
        b = nodeidx(fp, edge.nodeto)
        rc = edge.reducedcost
        fl = edge.flow
        return "$a=>$b(rc=$rc, flow=$fl)"
    end

    #println("Updating prices...")
    gamma = typemax(Int)

    i = fp.firstL
    while i !== nothing
        if i.inS

            i_idx = nodeidx(fp, i)

            #println("Edges to $(i_idx):")
            #printlist(i, :firstto, :nextto, printedge)
            #print("Active: ")
            #printlist(i, :firstactiveto, :nextactiveto, printedge)
            #print("Balanced: ")
            #printlist(i, :firstbalancedto, :nextbalancedto, printedge)
            #print("Inactive: ")
            #printlist(i, :firstinactiveto, :nextinactiveto, printedge)

            #println("Edges from $(i_idx):")
            #printlist(i, :firstfrom, :nextfrom, printedge)
            #print("Active: ")
            #printlist(i, :firstactivefrom, :nextactivefrom, printedge)
            #print("Balanced: ")
            #printlist(i, :firstbalancedfrom, :nextbalancedfrom, printedge)
            #print("Inactive: ")
            #printlist(i, :firstinactivefrom, :nextinactivefrom, printedge)

        end
        i = i.nextL
    end

    # TODO: Search with S or notS, based on set sizes (or count of edges in set?)
    i = fp.firstL
    while i !== nothing

        if i.inS

            # Adjust flows on balanced lines in to / out of S
            ij = i.firstbalancedfrom
            while ij !== nothing
                !ij.nodeto.inS && (ij.flow = ij.limit)
                ij = ij.nextbalancedfrom
            end

            ji = i.firstbalancedto
            while ji !== nothing
                !ji.nodefrom.inS && (ji.flow = 0)
                ji = ji.nextbalancedto
            end

            # Determine gamma from active / inactive lines in to / out of S
            ij = i.firstinactivefrom
            while ij !== nothing
                !ij.nodeto.inS && (gamma = min(gamma, ij.reducedcost))
                ij = ij.nextinactivefrom
            end

            ji = i.firstactiveto
            while ji !== nothing
                !ji.nodefrom.inS && (gamma = min(gamma, -ji.reducedcost))
                ji = ji.nextactiveto
            end

        end

        i = i.nextL

    end

    gamma === typemax(Int) && error("gamma === typemax(Int)")

    # Flows have changed, so recalculate imbalances
    # TODO: Do this more intelligently on the fly
    calculateimbalances!(fp)

    # Adjust node prices and edge reduced costs
    # TODO: Use S
    i = fp.firstL
    while i !== nothing

        if i.inS # TODO: Can do better, right now Ss are not contiguous!
                 # Need to add new Ls to opposite end of list from S selection

            # Price is increase by gamma
            i.price += gamma

            # Edges from i have reduced cost decreased by gamma
            ij = i.firstfrom
            while ij !== nothing

                j = ij.nodeto

                if !j.inS # S->S edges don't change

                    oldreducedcost = ij.reducedcost 
                    newreducedcost = oldreducedcost - gamma
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

                ij = ij.nextfrom

            end

            # Edges to i have reduced cost increased by gamma
            ji = i.firstto
            while ji !== nothing

                j = ji.nodefrom

                if !j.inS # S->S edges don't change

                    oldreducedcost = ji.reducedcost
                    newreducedcost = oldreducedcost + gamma
                    ji.reducedcost = newreducedcost

                    if oldreducedcost === 0 # ji moved from balanced -> inactive

                        # Remove edge from j's balancedfrom adjacency list
                        remove!(ji, :prevbalancedfrom, :nextbalancedfrom,
                                j, :firstbalancedfrom, :lastbalancedfrom)

                        # Remove edge from i's balancedto adjacency list
                        remove!(ji, :prevbalancedto, :nextbalancedto,
                                i, :firstbalancedto, :lastbalancedto)

                        # Add edge to j's inactivefrom adjacency list
                        addend!(ji, :previnactivefrom, :nextinactivefrom,
                                j, :firstinactivefrom, :lastinactivefrom)

                        # Add edge to i's inactiveto adjacency list
                        addend!(ji, :previnactiveto, :nextinactiveto,
                                i, :firstinactiveto, :lastinactiveto)

                    elseif newreducedcost === 0 # ji moved from active -> balanced

                        # Remove edge from j's activefrom adjacency list
                        remove!(ji, :prevactivefrom, :nextactivefrom,
                                j, :firstactivefrom, :lastactivefrom)

                        # Remove edge from i's activeto adjacency list
                        remove!(ji, :prevactiveto, :nextactiveto,
                                i, :firstactiveto, :lastactiveto)

                        # Add edge to j's balancedfrom adjacency list
                        addend!(ji, :prevbalancedfrom, :nextbalancedfrom,
                                j, :firstbalancedfrom, :lastbalancedfrom)

                        # Add edge to i's balancedto adjacency list
                        addend!(ji, :prevbalancedto, :nextbalancedto,
                                i, :firstbalancedto, :lastbalancedto)

                    end

                end

                ji = ji.nextto

            end

        end

        i = i.nextL

    end

    # Prices have changed, so recalculate ascent gradient
    # TODO: Do this more intelligently on the fly
    calculateascentgradient!(fp)

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
    # TODO: Use S instead
    i = fp.firstL
    while i !== nothing

        if i.inS

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

        end

        i = i.nextL

    end

end
