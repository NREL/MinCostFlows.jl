function solveflows!(fp::FlowProblem)

    while true

        # Find a potential starting point for an augmenting path
        augmentingpathstart = findfirst(x -> x > 0, fp.imbalances)

        # No starting point found, problem is either solved or infeasible
        augmentingpathstart == nothing && return fp

        # Reset the scan sets for the new iteration
        resetSL!(fp, augmentingpathstart)

        # Update either flows or shadow prices
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
function resetSL!(prob::FlowProblem, j::Int)

    for i in 1:prob.S
        S[i] && (S[i] = false)
        L[i] && (L[i] = false)
    end

    L[j] = true

    return prob

end

"""
If S is a subset of L, return an element of L - S after adding that element
to S. If S == L, return `nothing` and make no changes.
"""
function augmentS!(prob::FlowProblem)

    for i in 1:prob.n

        if prob.L[i] && !prob.S[i]
            prob.S[i] == true
            return i
        end

    end

    return nothing

end

"""
Checks whether the current set S would increase the
dual objective function if prices were updated
"""
function dualascendable(prob::FlowProblem)

    x = 0

    for i in 1:prob.n
        # TODO: Incorporate the other 4 terms in x
        if prob.S[i]
            x += prob.injections[i]
        end
    end

    return x > 0

end

"""
Add nodes to L that are adjacent to `i` and satisfy inclusion criteria,
returning one such added node that has a positive imbalance
(or `nothing` if none exist). Subset of Step 2 in Bertsekas.
"""
function augmentL!(prob::FlowProblem, i::Int)
    # TODO
    return nothing
end

"""
Updates the flows on the augmenting path (Step 3 in Bertsekas)
"""
function updateflows!(fp::FlowProblem, pathstart::Int, pathend::Int)
    # TODO
    return nothing
end

"""
Updates the shadow prices of the elements of S (Step 4 in Bertsekas)
"""
function updateprices!(fp::FlowProblem)
    # TODO
    return nothing
end
