struct LPFlowSolution
    flows::Vector{Int}
    prices::Vector{Float64}
end
 
"""
Generate a random flow problem on a network of size `(n, e)`,
with feasiblity fallback node.
"""
function randomproblem(n::Int, e::Int)

    # Create random connections / costs / limits / injections
    nodesfrom = rand(1:n, e)
    nodesto = nodesfrom .+ rand(1:(n-1), e)
    map!(i -> i > n ? i - n : i, nodesto, nodesto)
    flowlimits = rand(1:20, e)
    flowcosts = rand(1:5, e)
    injections = rand(-20:20, n)

    # Add fallback node
    push!(injections, -sum(injections))

    # Add edges to fallback node
    append!(nodesfrom, 1:n)
    append!(nodesto, fill(n+1, n))
    append!(flowlimits, fill(9999, n))
    append!(flowcosts, fill(0, n))

    # Add fallback edges from fallback node
    append!(nodesfrom, fill(n+1, n))
    append!(nodesto, 1:n)
    append!(flowlimits, fill(9999, n))
    append!(flowcosts, fill(9999, n)) 

    return FlowProblem(nodesfrom, nodesto, flowlimits, flowcosts, injections)

end

"""
Formulate and solve an LP based on a FlowProblem.
"""
function MathProgBase.linprog(fp::FlowProblem)
    
    A = buildAmatrix(fp)

    #display(full(A)); println()
    #display(fp.injections); println()
    #display(fp.limits); println()
    #display([full(A) .-fp.injections])

    # Solve the LP
    # println("Net injection: ", sum(fp.injections))
    # display(fp.injections); println()
    result = linprog(Vector{Float64}(fp.costs), A, '=', Vector{Float64}(.-fp.injections),
                     0., Vector{Float64}(fp.limits), ClpSolver())
    (result.status != :Optimal) && error("LP did not solve")
    return LPFlowSolution(result.sol, result.attrs[:lambda])

end

function buildAmatrix(fp::FlowProblem)

    n_edges = length(fp.nodesfrom)
    n_nodes = fp.nodes

    Is = Vector{Float64}(undef, 2*n_edges)
    Js = similar(Is)
    Vs = similar(Is)

    # Construct the sparse A matrix
    for ij in 1:n_edges

        i = fp.nodesfrom[ij]
        j = fp.nodesto[ij]

        Is[2ij-1] = i
        Js[2ij-1] = ij
        Vs[2ij-1] = -1

        Is[2ij] = j
        Js[2ij] = ij
        Vs[2ij] = 1
        
    end

    return sparse(Is, Js, Vs)

end
