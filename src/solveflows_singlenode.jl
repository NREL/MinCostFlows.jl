"""
Returns true if the process updated flows or prices, and false otherwise.
"""
function singlenodeupdate!(fp::FlowProblem, i::Node)

    return false
    statechanged = false

    while true

        # Determine B+ and B-
        # Step 1: Compare i.imbalance with B+ and B- room
        x = 0
        for edge in Bplus
            x += edge.limit
            x -= edge.flow
        end
        for edge in Bminus
            x += edge.flow
        end

        if i.imbalance >= x

            # Step 4: Max/min out flows and increase i.price
            updateprice!(fp, i)
            statechanged = true
            i.imbalance > 0 ? continue : break

        elseif i.imbalance > 0 && length(Bplus) > 0

            # Step 2: Outgoing arc flow adjustment
            updateoutgoingflows!(fp, i)
            statechanged = true
            continue

        elseif i.imbalance > 0 && length(Bminus) > 0

            # Step 3: Incoming arc flow adjustment
            updateincomingflows!(fp, i)
            statechanged = true
            continue

        else # No nodes in B, or i.imbalance == 0
            break

        end

    end

    return statechanged

end

function updateoutgoingflows!(fp::FlowProblem, i::Node)

end

function updateincomingflows!(fp::FlowProblem, i::Node)

end

function updateprice!(fp::FlowProblem, i::Node)

end
