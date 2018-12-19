abstract type AbstractNode end

mutable struct Edge{T<:AbstractNode}
    nodefrom::T
    nodeto::T
    nextfrom::Union{Edge{T},Nothing}
    nextto::Union{Edge{T},Nothing}
    previntoS::Union{Edge{T},Nothing}
    nextintoS::Union{Edge{T},Nothing}
    prevoutofS::Union{Edge{T},Nothing}
    nextoutofS::Union{Edge{T},Nothing}
    limit::Int
    cost::Int
    reducedcost::Int
    flow::Int
    forward::Bool
    backward::Bool # TODO: Switch to linked list and eliminate

    Edge{}(from::T, to::T, limit::Int, cost::Int) where {T <: AbstractNode} =
        new{T}(from, to, nothing, nothing,
               nothing, nothing, nothing, nothing, limit, cost, cost, 0, false, false)
end

mutable struct Node <: AbstractNode
    firstfrom::Union{Edge{Node},Nothing}
    firstto::Union{Edge{Node},Nothing}
    label::Union{Edge{Node},Nothing}
    nextL::Union{Node,Nothing}
    injection::Int
    price::Int
    imbalance::Int
    inS::Bool
    inL::Bool

    Node(injection::Int) = new(nothing, nothing, nothing, nothing,
                               injection, 0, injection, false, false)
end

mutable struct FlowProblem
   nodes::Vector{Node}
   edges::Vector{Edge{Node}}
   firstL::Union{Node,Nothing}
   firstintoS::Union{Edge{Node},Nothing}
   firstoutofS::Union{Edge{Node},Nothing}
end

function FlowProblem(nodesfrom::Vector{Int}, nodesto::Vector{Int},
                     limits::Vector{Int}, costs::Vector{Int},
                     injections::Vector{Int})

    n = length(injections)
    e = length(nodesfrom)

    @assert length(nodesto) == e
    @assert length(limits) == e
    @assert length(costs) == e
    @assert sum(injections) == 0

    nodes = Node.(injections)
    edges = Edge.(nodes[nodesfrom], nodes[nodesto], limits, costs)

    for (e, edge) in enumerate(edges)

        nodefrom = edge.nodefrom
        nodeto = edge.nodeto

        nodefrom.firstfrom === nothing && (nodefrom.firstfrom = edge)
        nodeto.firstto === nothing && (nodeto.firstto = edge)

        edge.nextfrom = next(x -> x.nodefrom === nodefrom, edges, e+1)
        edge.nextto = next(x -> x.nodeto === nodeto, edges, e+1)

    end

    return FlowProblem(nodes, edges, nothing, nothing, nothing)

end

flows(fp::FlowProblem) = getproperty.(fp.edges, :flow)
costs(fp::FlowProblem) = getproperty.(fp.edges, :cost)
limits(fp::FlowProblem) = getproperty.(fp.edges, :limit)
injections(fp::FlowProblem) = getproperty.(fp.nodes, :injection)
prices(fp::FlowProblem) = getproperty.(fp.nodes, :price)

function next(f::Function, v::Vector, start::Int)

    for i in start:length(v)
        x = v[i]
        f(x) && return x
    end

    return nothing

end
