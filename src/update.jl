"""
Change injections on `targetnode` to `newinj` for a FlowProblem `fp`,
automatically adjusting `slacknode` to compensate.
"""
function updateinjection!(targetnode::Node, slacknode::Node, newinj::Int)

    change = newinj - targetnode.injection

    targetnode.injection += change
    targetnode.imbalance += change
    slacknode.injection -= change
    slacknode.imbalance -= change

    return

end

"""
Change line limit on `edge` to `newlimit` for a FlowProblem `fp`.
"""
function updateflowlimit!(edge::Edge{Node}, newlimit::Int)

    newlimit < 0 && error("New limit must be non-negative, received $newlimit")

    if (newlimit < edge.flow) || (edge.reducedcost < 0)
        change = edge.flow - newlimit
        edge.flow = newlimit
        edge.nodefrom.imbalance += change
        edge.nodeto.imbalance -= change
    end

    edge.limit = newlimit

    return

end
