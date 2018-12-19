"""
Implements the Relaxation dual-ascent method for solving
min-cost flow problems, as described in Bertsekas (1998)
"""
function solveflows!(fp::FlowProblem)

    elist = edgelist(fp)

    recalculateimbalances!(fp)

    while true

        println(elist)
        println("Costs: ", costs(fp))
        println("Limits: ", limits(fp))
        println("Flows: ", flows(fp))
        println("Injections: ", injections(fp))
        println("Prices: ", prices(fp))
        @assert complementarityslackness(fp)
        #println("Starting major iteration")

        # Find a potential starting node for an augmenting path
        augmentingpathstart = firstpositiveimbalance(fp)

        # No starting node found, problem is either solved or infeasible
        augmentingpathstart === nothing && return fp

        # Reset the scan sets for this iteration
        resetSL!(fp, augmentingpathstart)
        lengthSinout(fp)
        #showSL(fp)

        # Main iteration: update either flows or shadow prices
        update!(fp, augmentingpathstart)

    end

end

function update!(fp::FlowProblem, augmentingpathstart::Node)

    while true

        #println("Starting minor iteration")

        # Look for a candidate node i to scan and add to S
        i = augmentS!(fp)
        lengthSinout(fp)
        #showSL(fp)

        # Update prices if it will improve the dual solution
        # or if there are no nodes left to scan
        (i == nothing || dualascendable(fp)) && return updateprices!(fp)

        # Label neighbour nodes of i
        augmentingpathend = augmentL!(fp, i)
        #showSL(fp)

        # Didn't find an augmenting path, try adding a different node
        augmentingpathend == nothing && continue

        # Found an augmenting path, augment flows accordingly
        return updateflows!(fp, augmentingpathstart, augmentingpathend)

    end

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

    # Reset into/outof edge lists
    fp.firstintoS = nothing
    fp.firstoutofS = nothing

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

            ij = i.firstfrom
            while ij !== nothing
                if ij.nodeto.inS # i and j are both in S now
                    # Extract ij from the intoS linked list
                    if ij.previntoS === nothing # ij is the first element of the list
                        fp.firstintoS = ij.nextintoS
                    else # ij is not the first element of the list
                        ij.previntoS.nextintoS = ij.nextintoS
                    end
                else # ij leads out of S
                    # Add ij to the outofS linked list
                    if fp.firstoutofS === nothing # ij is the only element
                        ij.nextoutofS = nothing
                    else # ij goes to the beginning of an existing list
                        fp.firstoutofS.prevoutofS = ij
                        ij.nextoutofS = fp.firstoutofS
                    end
                    ij.prevoutofS = nothing
                    fp.firstoutofS = ij
                end
                ij = ij.nextfrom
            end

            ji = i.firstto
            while ji !== nothing
                if ji.nodefrom.inS # j and i are both in S now
                    # Extract ji from the outofS linked list
                    if ji.prevoutofS === nothing # ji is the first element of the list
                        fp.firstoutofS = ji.nextoutofS
                    else # ij is not the first element of the list
                        ji.prevoutofS.nextoutofS = ji.nextoutofS
                    end 
                else # ji leads in to S
                    # Add ji to the intoS linked list
                    if fp.firstintoS === nothing # ji is the only element
                        ji.nextintoS = nothing
                    else # ji goes to the beginning of an existing list
                        fp.firstintoS.prevoutofS = ji
                        ji.nextintoS = fp.firstintoS
                    end
                    ji.previntoS = nothing
                    fp.firstintoS = ji
                end
                ji = ji.nextto
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

    #println("Calculating ascent gradient...")
    x = 0

    node = fp.firstL
    while node !== nothing
        node.inS && (x += node.injection) #TODO: Fix S ordering in L list to avoid this
        node = node.nextL
    end

    ij = fp.firstintoS
    while ij !== nothing
        # ij leads in to S and is active
        (ij.reducedcost < 0) && (x += ij.limit)
        ij = ij.nextintoS
    end

    ij = fp.firstoutofS
    while ij !== nothing
        # ij leads out of S and is balanced or active
        (ij.reducedcost <= 0) && (x -= ij.limit)
        ij = ij.nextoutofS
    end

    #println("Gradient: ", x)
    return x > 0

end

"""
Add nodes to L that are adjacent to `i` and satisfy inclusion criteria,
returning one such added node that has a positive imbalance
(or `nothing` if none exist). Subset of Step 2 in Bertsekas.
"""
function augmentL!(fp::FlowProblem, i::N)::Union{N,Nothing} where {N<:Node}

    #println("Augmenting L...")
    negativeimbalance = nothing

    # Iterate through edges ji connecting j to i
    ji = i.firstto
    while ji !== nothing

        j = ji.nodefrom # source node of edge ji

        if (!j.inL && # j isn't already in L #TODO: Linked list again?
            ji.reducedcost === 0 && # ji is balanced
            ji.flow > 0) # flow exceeds lower bound

            # Add j to L
            j.inL = true
            j.nextL = fp.firstL
            fp.firstL = j

            # Define label for j
            j.label = ji

            # Check for negative imbalance and save for returning
            j.imbalance < 0 && (negativeimbalance = j)

        end

        ji = ji.nextto # Get next edge

    end

    # Iterate through edges ij connecting i to j
    ij = i.firstfrom
    while ij !== nothing

        j = ij.nodeto # sink node of edge ij

        if (!j.inL && # j isn't already in L
            ij.reducedcost === 0 && # ij is balanced
            (ij.flow < ij.limit)) # ij flow is less than limit

            # Add j to L
            j.inL = true
            j.nextL = fp.firstL
            fp.firstL = j

            # Define label for j
            j.label = ij

            # Check for negative imbalance and save for returning
            j.imbalance < 0 && (negativeimbalance = j)

        end

        ij = ij.nextfrom # Get next edge

    end

    return negativeimbalance

end

"""
Updates the flows on the augmenting path from `startnode` to `endnode`
(Step 3 in Bertsekas)
"""
function updateflows!(fp::FlowProblem, startnode::Node, endnode::Node)

    #println("Updating flows...")
    currentnode = endnode
    delta = min(startnode.imbalance, -endnode.imbalance)
    for edge in fp.edges # TODO: Eliminate with linked list?
        edge.forward = false
        edge.backward = false
    end

    # Traverse the labels backwards to startnode
    while currentnode !== startnode

        # Determine edge direction and move to previous node
        ab = currentnode.label # Get edge
        a = ab.nodefrom # Get edge to/from nodes
        b = ab.nodeto

        if b === currentnode # ab is a forward edge in the path, move to node a
            ab.forward = true
            currentnode = a
            delta = min(delta, ab.limit - ab.flow)
        else # ab is a backwards edge in the path, move to node b
            ab.backward = true
            currentnode = b
            delta = min(delta, ab.flow)
        end

    end

    # Adjust edge flows as appropriate
    for edge in fp.edges # TODO: Switch to linked list
        edge.forward && (edge.flow += delta)
        edge.backward && (edge.flow -= delta)
    end

    recalculateimbalances!(fp) #TODO: Do this on the fly instead

end

"""
Updates the shadow prices (and potentially flows) of the elements of S
(Step 4 in Bertsekas)
"""
function updateprices!(fp::FlowProblem)

    #println("Updating prices...")
    gamma = typemax(Int)

    ij = fp.firstoutofS
    while ij !== nothing
      if ij.reducedcost === 0
          ij.flow = ij.limit # Adjust flow
      elseif ij.flow < ij.limit
          gamma = min(gamma, ij.reducedcost)
      end
      ij = ij.nextoutofS
    end

    ij = fp.firstintoS
    while ij !== nothing
        if ij.reducedcost === 0
            ij.flow = 0
        elseif ij.flow > 0
            gamma = min(gamma, -ij.reducedcost)
        end
      ij = ij.nextintoS
    end

    # Adjust node prices and edge reduced costs
    i = fp.firstL
    while i !== nothing

        if i.inS # TODO: Can do better, right now Ss are not contiguous!
                 # Need to add new Ls to opposite end of list from S selection

            # Price is increase by gamma
            i.price += gamma

            # Edges from i have reduced cost decreased by gamma
            ij = i.firstfrom
            while ij !== nothing
                ij.reducedcost -= gamma
                ij = ij.nextfrom
            end

            # Edges to i have reduced cost increased by gamma
            ji = i.firstto
            while ji !== nothing
                ji.reducedcost += gamma
                ji = ji.nextto
            end

        end

        i = i.nextL

    end

    recalculateimbalances!(fp)

end

# TODO: Should probably recalculate these on the fly instead (whenever flow is updated)
function recalculateimbalances!(fp::FlowProblem)

    for node in fp.nodes
        node.imbalance = node.injection
    end

    for edge in fp.edges
        edge.nodefrom.imbalance -= edge.flow
        edge.nodeto.imbalance += edge.flow
    end

end
