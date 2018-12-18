"""
Implements the Relaxation dual-ascent method for solving
min-cost flow problems, as described in Bertsekas (1998)
"""
function solveflows!(fp::FlowProblem)

    recalculateimbalances!(fp)

    while true

        # Find a potential starting point for an augmenting path
        augmentingpathstart = findfirst(x -> x > 0, fp.imbalances)

        # No starting point found, problem is either solved or infeasible
        augmentingpathstart == nothing && return fp

        # Reset the scan sets for this iteration
        resetSL!(fp, augmentingpathstart)

        # Main iteration: update either flows or shadow prices
        update!(fp, augmentingpathstart)

    end

end

function update!(fp::FlowProblem, augmentingpathstart::Int)

    while true

        # Look for a candidate node i to scan and add to S
        i = augmentS!(fp)

        # Update prices if it will improve the dual solution
        # or if there are no nodes left to scan
        (i == nothing || dualascendable(fp)) && return updateprices!(fp)

        # Label neighbour nodes of i
        augmentingpathend = augmentL!(fp, i)

        # Didn't find an augmenting path, try adding a different node
        augmentingpathend == nothing && continue

        # Found an augmenting path, augment flows accordingly
        return updateflows!(fp, augmentingpathstart, augmentingpathend)

    end

end

"""
Empties the set S and reduces the set L to a single element, `j`
"""
function resetSL!(fp::FlowProblem, j::Int)

    fill!(fp.S, false)
    fill!(fp.transS, 0)

    fill!(fp.L, false)
    # TODO: Consider using transL?
    fp.L[j] = true

    return nothing

end

"""
If S is a subset of L, return an element of L - S after adding that element
to S. If S == L, return `nothing` and make no changes.
"""
function augmentS!(fp::FlowProblem)

    for i in 1:fp.nodes

        if fp.L[i] && !fp.S[i]

            # i joins S
            fp.S[i] = true

            ij = fp.firstfrom[i]
            while ij > 0
                if fp.S[fp.nodesto[ij]] # i and j are both in S now
                    fp.transS[ij] = 0
                else # ij leads out of S
                    fp.transS[ij] = -1
                end
                ij = fp.nextfrom[ij]
            end

            ji = fp.firstto[i]
            while ji > 0
                if fp.S[fp.nodesfrom[ji]] # i and j are both in S now
                    fp.transS[ji] = 0
                else # ji leads in to S
                    fp.transS[ji] = 1
                end
                ji = fp.nextto[ji]
            end

            return i

        end

    end

    return nothing

end

"""
Checks whether the current set S would increase the
dual objective function if prices were updated
"""
function dualascendable(fp::FlowProblem)

    x = 0

    @inbounds for i in 1:fp.nodes
        fp.S[i] && (x += fp.injections[i])
    end

    @inbounds for ij in 1:length(fp.nodesfrom)

        (fp.transS[ij] === 0) && continue

        # ij leads in to S and is active
        if (fp.transS[ij] === 1) && (fp.reducedcosts[ij] < 0) # Spends a lot of time promoting?
            x += fp.limits[ij]

        # ij leads out of S and is balanced or active
        elseif fp.reducedcosts[ij] <= 0 # Spends a lot of time promoting?
            x -= fp.limits[ij]

        end

    end

    return x > 0

end

"""
Add nodes to L that are adjacent to `i` and satisfy inclusion criteria,
returning one such added node that has a positive imbalance
(or `nothing` if none exist). Subset of Step 2 in Bertsekas.
"""
function augmentL!(fp::FlowProblem, i::Int)

    negativeimbalance = 0

    # Iterate through edges ji connecting j to i
    ji = fp.firstto[i]
    while ji != 0

        j = fp.nodesfrom[ji] # source node of edge ji

        if (!fp.L[j] && # j isn't already in L
           isbalanced(fp, j, ji, i) && # ji is balanced
           fp.flows[ji] > 0) # flow exceeds lower bound

            # Add j to L
            fp.L[j] = true

            # Define label for j
            fp.labels[j] = ji

            # Check for negative imbalance and save for returning
            fp.imbalances[j] < 0 && (negativeimbalance = j)

        end

        ji = fp.nextto[ji] # Get next edge

    end

    # Iterate through edges ij connecting i to j
    ij = fp.firstfrom[i]
    while ij != 0

        j = fp.nodesto[ij] # sink node of edge ij

        if (!fp.L[j] && # j isn't already in L
           isbalanced(fp, i, ij, j) && # ij is balanced
           (fp.flows[ij] < fp.limits[ij])) # ij flow is less than limit

            # Add j to L
            fp.L[j] = true

            # Define label for j
            fp.labels[j] = ij

            # Check for negative imbalance and save for returning
            fp.imbalances[j] < 0 && (negativeimbalance = j)

        end

        ij = fp.nextfrom[ij] # Get next edge

    end

    return negativeimbalance > 0 ? negativeimbalance : nothing

end

"""
Updates the flows on the augmenting path from `startnode` to `endnode`
(Step 3 in Bertsekas)
"""
function updateflows!(fp::FlowProblem, startnode::Int, endnode::Int)

    currentnode = endnode
    delta = min(fp.imbalances[startnode], -fp.imbalances[endnode])
    fill!(fp.forwardedges, false)
    fill!(fp.backwardedges, false)

    # Traverse the labels backwards to startnode
    while currentnode != startnode

        # Determine edge direction and move to previous node
        ab = fp.labels[currentnode] # Get edge number
        a = fp.nodesfrom[ab] # Get edge to/from nodes
        b = fp.nodesto[ab]

        if b == currentnode # ab is a forward edge in the path, move to node a
            fp.forwardedges[ab] = true
            currentnode = a
            deltacandidate = fp.limits[ab] - fp.flows[ab]
        else # ab is a backwards edge in the path, move to node b
            fp.backwardedges[ab] = true
            currentnode = b
            deltacandidate = fp.flows[ab]
        end

        # Decrease delta if nessecary
        delta = min(delta, deltacandidate)

    end

    # Adjust edge flows as appropriate
    for i in 1:length(fp.flows)
        fp.forwardedges[i] && (fp.flows[i] += delta)
        fp.backwardedges[i] && (fp.flows[i] -= delta)
    end
    recalculateimbalances!(fp)

    return nothing

end

"""
Updates the shadow prices (and potentially flows) of the elements of S
(Step 4 in Bertsekas)
"""
function updateprices!(fp::FlowProblem)

    gamma = typemax(Int)

    # TODO: Convert to use transS
    for i in 1:length(fp.S)

        fp.S[i] || continue # Only proceed if i is in S

        # Cycle through edges ij, determine if edge j is outside of S
        ij = fp.firstfrom[i]
        while ij != 0

            j = fp.nodesto[ij]

            if !fp.S[j] # j is outside of S
                if isbalanced(fp, i, ij, j)
                    fp.flows[ij] = fp.limits[ij] # Adjust flow
                elseif fp.flows[ij] < fp.limits[ij]
                    gamma = min(gamma, fp.shadowprices[j] + fp.costs[ij] - fp.shadowprices[i])
                end
            end

            ij = fp.nextfrom[ij]

        end

        # Cycle through edges ji, determine if edge j is outside of S
        ji = fp.firstto[i]
        while ji != 0

            j = fp.nodesfrom[ji]

            if !fp.S[j] # j is outside of S
                if isbalanced(fp, j, ji, i)
                    fp.flows[ji] = 0 # Adjust flow
                elseif fp.flows[ji] > 0
                    gamma = min(gamma, fp.shadowprices[j] - fp.costs[ji] - fp.shadowprices[i])
                end
            end

            ji = fp.nextto[ji]

        end

    end

    # Adjust node prices and edge reduced costs
    for i in 1:length(fp.S)
        if fp.S[i]

            # Price is increase by gamma
            fp.shadowprices[i] += gamma

            # Edges from i have reduced cost decreased by gamma
            ij = fp.firstfrom[i]
            while ij > 0
                fp.reducedcosts[ij] -= gamma
                ij = fp.nextfrom[ij]
            end

            # Edges to i have reduced cost increased by gamma
            ji = fp.firstto[i]
            while ji > 0
                fp.reducedcosts[ji] += gamma
                ji = fp.nextto[ji]
            end

        end
    end

    recalculateimbalances!(fp)
    return nothing

end

isinactive(fp::FlowProblem, i::Int, ij::Int, j::Int) =
    fp.shadowprices[i] < fp.costs[ij] + fp.shadowprices[j]

isbalanced(fp::FlowProblem, i::Int, ij::Int, j::Int) =
    fp.shadowprices[i] == fp.costs[ij] + fp.shadowprices[j]

isactive(fp::FlowProblem, i::Int, ij::Int, j::Int) =
    fp.shadowprices[i] > fp.costs[ij] + fp.shadowprices[j]

# TODO: Should probably recalculate these on the fly (whenever flow is updated)
function recalculateimbalances!(fp::FlowProblem)

    fp.imbalances .= fp.injections

    for ij in 1:length(fp.nodesfrom)

        flow = fp.flows[ij]

        i = fp.nodesfrom[ij]
        fp.imbalances[i] -= flow

        j = fp.nodesto[ij]
        fp.imbalances[j] += flow

    end

end
