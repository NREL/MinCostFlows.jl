function multinodeupdate!(fp::FlowProblem, augmentingpathstart::Node)

    iters = 0

    # Reset the scan sets for this iteration
    resetSL!(fp, augmentingpathstart)
    #showSL(fp)

    while true

        iters += 1
        #println("Starting minor iteration")

        # Look for a candidate node i to scan and add to S
        i = augmentS!(fp) # ~ 1/3 of time here
        #showSL(fp)

        # Update prices if it will improve the dual solution
        # or if there are no nodes left to scan
        if i === nothing || dualascendable(fp)
            updateprices!(fp)
            break
        end

        # Label neighbour nodes of i
        augmentingpathend = augmentL!(fp, i) # ~ 2/3 of time here
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
Empties the set S and reduces the set L to a single element, `j`
"""
function resetSL!(fp::FlowProblem, j::Node)

    #println("Resetting SL...")

    node = fp.firstS
    while node !== nothing
        node.inL = false
        node.inS = false
        node = node.nextS
    end

    node = fp.firstLnotS
    while node !== nothing
        node.inL = false
        node = node.nextLnotS
    end

    # Reset S to {}
    fp.firstS = nothing
    fp.lastS = nothing
    fp.ascentgradient = 0

    # Reset L to {j}
    fp.firstLnotS = j
    j.prevLnotS = nothing
    j.nextLnotS = nothing
    fp.lastLnotS = j
    j.inL = true

    return nothing

end

"""
If S is a subset of L, return an element of L - S after adding that element
to S. If S == L, return `nothing` and make no changes.
"""
function augmentS!(fp::FlowProblem)

    #println("Augmenting S...")
    i = fp.firstLnotS
    i === nothing && return nothing

    # i joins S and leaves LnotS
    i.inS = true
    addend!(i, :prevS, :nextS, fp, :firstS, :lastS)
    remove!(i, :prevLnotS, :nextLnotS, fp, :firstLnotS, :lastLnotS)

    # Update the ascent gradient with i included in S
    ij = i.firstbalancedfrom
    while ij !== nothing
        # TODO: Could factor out ij.flow here (save an Int allocation?)
        if ij.nodeto.inS
            fp.ascentgradient += ij.flow
        else
            fp.ascentgradient -= (ij.limit - ij.flow)
        end
        ij = ij.nextbalancedfrom
    end

    ji = i.firstbalancedto
    while ji !== nothing
        # TODO: Could factor out ji.flow here (save an Int allocation?)
        if ji.nodefrom.inS
            fp.ascentgradient += (ji.limit - ji.flow)
        else
            fp.ascentgradient -= ji.flow
        end
        ji = ji.nextbalancedto
    end

    return i

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

            # Add j to LnotS
            j.inL = true
            addend!(j, :prevLnotS, :nextLnotS, fp, :firstLnotS, :lastLnotS)

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

            # Add j to LnotS
            j.inL = true
            addend!(j, :prevLnotS, :nextLnotS, fp, :firstLnotS, :lastLnotS)

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

    #function printedge(edge::Edge)
    #    a = nodeidx(fp, edge.nodefrom)
    #    b = nodeidx(fp, edge.nodeto)
    #    rc = edge.reducedcost
    #    fl = edge.flow
    #    return "$a=>$b(rc=$rc, flow=$fl)"
    #end

    #println("Updating prices...")
    gamma = typemax(Int)

    #i = fp.firstS
    #while i !== nothing

        #i_idx = nodeidx(fp, i)

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

        #i = i.nextS

    #end

    # TODO: Search with S or notS, based on set sizes
    i = fp.firstS
    while i !== nothing


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

        i = i.nextS

    end

    gamma === typemax(Int) && error("gamma === typemax(Int)")

    # Flows have changed, so recalculate imbalances
    # TODO: Do this more intelligently on the fly
    calculateimbalances!(fp)

    # Adjust node prices and edge reduced costs
    i = fp.firstS
    while i !== nothing

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

        i = i.nextS

    end

    # Prices have changed, so recalculate ascent gradient
    # TODO: Do this more intelligently on the fly
    calculateascentgradient!(fp)

end
