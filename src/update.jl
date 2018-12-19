"""
Change injections on `targetnode` to `newinj` for a FlowProblem `fp`,
automatically adjusting `slacknode` to compensate.
"""
function updateinjection!(targetnode::Node, slacknode::Node, newinj::Int)

    change = newinj - targetnode.injection

    targetnode.injection += change
    slacknode.injection -= change

    return

end

"""
Change line limit on `edge` to `value` for a FlowProblem `fp`.
"""
function updateflowlimit!(edge::Edge{Node}, newlimit::Int)

    if (newlimit < edge.flow) || (edge.reducedcost < 0)
        edge.flow = newlimit
    end

    edge.limit = newlimit

    return

end
