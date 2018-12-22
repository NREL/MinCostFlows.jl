# Singly-linked lists

function addstart!(elem, next::Symbol, context, first::Symbol)

    oldfirstelem = getproperty(context, first)

    if oldfirstelem === nothing # list is empty

        # new element is last as well as first
        setproperty!(elem, next, nothing)

    else # other elements exist

        # add elem before current first element
        setproperty!(elem, next, oldfirstelem)

    end

    # elem is now first element
    setproperty!(context, first, elem)

end

function remove!(context, first::Symbol, next::Symbol)

    oldfirstelem = getproperty(context, first)

    if oldfirstelem !== nothing # nothing to do for an empty list
        newfirstelem = getproperty(oldfirstelem, next)
        setproperty!(context, first, newfirstelem)
    end

    return

end

# Doubly-linked lists

function addstart!(elem, prev::Symbol, next::Symbol, context, first::Symbol, last::Symbol)

    oldfirstelem = getproperty(context, first)

    if oldfirstelem === nothing # list is empty

        # new element is last as well as first
        setproperty!(elem, next, nothing)
        setproperty!(context, last, elem)

    else # other elements exist

	# add elem before current first element
        setproperty!(elem, next, oldfirstelem)
        setproperty!(oldfirstelem, prev, elem)

    end

    # elem is now first element
    setproperty!(context, first, elem)
    setproperty!(elem, prev, nothing)

end

function addend!(elem, prev::Symbol, next::Symbol, context, first::Symbol, last::Symbol)

    oldlastelem = getproperty(context, last)

    if oldlastelem === nothing # list is empty

        # new element is first as well as last
        setproperty!(elem, prev, nothing)
        setproperty!(context, first, elem)

    else # other elements exist

        # add elem after current last element
        setproperty!(elem, prev, oldlastelem)
        setproperty!(oldlastelem, next, elem)

    end

    # elem is now last element
    setproperty!(context, last, elem)
    setproperty!(elem, next, nothing)

end

# e.g. remove!(ij, :prevbalancedto, :nextbalancedto, i, :firstbalancedto, :lastbalancedto)
function remove!(elem, prev::Symbol, next::Symbol, context, first::Symbol, last::Symbol)

    prevelem = getproperty(elem, prev)
    nextelem = getproperty(elem, next)

    if prevelem === nothing # elem is at start of list
        setproperty!(context, first, nextelem)
    else
        setproperty!(prevelem, next, nextelem)
    end

    if nextelem === nothing # elem is at end of list
        setproperty!(context, last, prevelem)
    else
        setproperty!(nextelem, prev, prevelem)
    end

    return

end
