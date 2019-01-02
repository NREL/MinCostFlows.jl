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
Change line limit on `edge` to `value` for a FlowProblem `fp`.
"""
function updateflowlimit!(edge::Edge{Node}, newlimit::Int)

    if (newlimit < edge.flow) || (edge.reducedcost < 0)
        change = edge.flow - newlimit
        edge.flow = newlimit
        edge.nodefrom.imbalance += change
        edge.nodeto.imbalance -= change
    end

    edge.limit = newlimit

    return

end
