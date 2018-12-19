defzero(::Nothing) = 0
defzero(x::Int) = x

function showSL(fp::FlowProblem)

    S = Int[]
    L = Int[]
    print("L: [")

    node = fp.firstL
    while node !== nothing
        i = findfirst(n -> n === node, fp.nodes)
        print(i, ", ")
        i in L && error("Loop in L!!")
        push!(L, i)
        node.inS && push!(S, i)
        node = node.nextL
    end

    println("]")
    println("S: ", S)

end
