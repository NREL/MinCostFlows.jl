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

    # Eliminate duplicate edges
    fromtos = Tuple{Int,Int}[]
    fromto_idxs = Int[]
    for i in 1:e
        fromto = (nodesfrom[i], nodesto[i])
        fromto in fromtos && continue
        push!(fromtos, fromto)
        push!(fromto_idxs, i)
    end

    e = length(fromto_idxs)
    nodesfrom = nodesfrom[fromto_idxs]
    nodesto = nodesto[fromto_idxs]

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

    result = linprog(Vector{Float64}(costs(fp)), A, '=', Vector{Float64}(.-injections(fp)),
                     0., Vector{Float64}(limits(fp)), ClpSolver())

    (result.status != :Optimal) && error("LP did not solve")

    return LPFlowSolution(result.sol, result.attrs[:lambda])

end

function buildAmatrix(fp::FlowProblem)

    n_edges = length(fp.edges)
    n_nodes = length(fp.nodes)

    Is = Vector{Int}(undef, 2*n_edges)
    Js = Vector{Int}(undef, 2*n_edges)
    Vs = Vector{Float64}(undef, 2*n_edges)

    # Construct the sparse A matrix
    for (ij, edge) in enumerate(fp.edges)

        i = edge.nodefrom
        j = edge.nodeto

        Is[2ij-1] = findfirst(x -> x === i, fp.nodes)
        Js[2ij-1] = ij
        Vs[2ij-1] = -1

        Is[2ij] = findfirst(x -> x === j, fp.nodes)
        Js[2ij] = ij
        Vs[2ij] = 1
        
    end

    return sparse(Is, Js, Vs)

end
