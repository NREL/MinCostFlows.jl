"""
Change injections on `targetnode` to `newinj` for a FlowProblem `fp`,
automatically adjusting `slacknode` to compensate.
"""
function updateinjection!(fp::FlowProblem, newinj::Int, targetnode::Int, slacknode::Int)
    change = newinj - fp.injections[targetnode]
    fp.injections[targetnode] += change
    fp.injections[slacknode] -= change
    return
end

"""
Change line limit on `edge` to `value` for a FlowProblem `fp`.
"""
function updateflowlimit!(fp::FlowProblem, newlimit::Int, edge::Int)

    if (newlimit < fp.flows[edge]) ||
        isactive(fp, fp.nodesfrom[edge], edge, fp.nodesto[edge])
        fp.flows[edge] = newlimit
    end

    fp.limits[edge] = newlimit

    return

end
