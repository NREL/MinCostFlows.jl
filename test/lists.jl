using InteractiveUtils

@testset "Link Operations" begin

    @testset "Singly-Linked Lists" begin

        ctx = Context()
        l1 = Link(1)
        l2 = Link(2)

        #function f(l::Link, c::Context)
        #    MinCostFlow.@addstart!(l, :nxt, c, :fst)
        #    MinCostFlow.@remove!(c, :fst, :nxt)
        #    return c.fst
        #end
        #@code_warntype f(l2, ctx)
        #println(@macroexpand(MinCostFlow.@addstart!(l2, :nxt, ctx, :fst)))

        @test ctx.fst === nothing
        
        MinCostFlow.@addstart!(l2, :nxt, ctx, :fst)
        MinCostFlow.@addstart!(l1, :nxt, ctx, :fst)

        @test ctx.fst === l1
        @test l1.nxt === l2
        @test l2.nxt === nothing

        @macroexpand MinCostFlow.@remove!(ctx, :fst, :nxt)
        MinCostFlow.@remove!(ctx, :fst, :nxt)
        @test ctx.fst === l2
        @test l2.nxt === nothing

    end

    @testset "Doubly-Linked Lists" begin

        ctx = Context()
        l1 = Link(1)
        l2 = Link(2)
        l3 = Link(3)
        l4 = Link(4)
        l5 = Link(5)

        #function f(la::Link, lb::Link, c::Context)
        #    MinCostFlow.@addstart!(la, :prv, :nxt, c, :fst, :lst)
        #    MinCostFlow.@addend!(lb, :prv, :nxt, c, :fst, :lst)
        #    MinCostFlow.@remove!(la, :prv, :nxt, c, :fst, :lst)
        #    return c.fst
        #end
        #@code_warntype f(l1, l2, ctx)

        @test ctx.fst === nothing
        @test ctx.lst === nothing

        # Build list
        MinCostFlow.@addstart!(l3, :prv, :nxt, ctx, :fst, :lst)
        MinCostFlow.@addend!(l4, :prv, :nxt, ctx, :fst, :lst)
        MinCostFlow.@addstart!(l2, :prv, :nxt, ctx, :fst, :lst)

        @test firstfromlast(ctx) === ctx.fst
        @test lastfromfirst(ctx) === ctx.lst

        @test ctx.fst === l2
        @test l2.prv === nothing
        @test l2.nxt === l3
        @test l3.prv === l2
        @test l3.nxt === l4
        @test l4.prv === l3
        @test l4.nxt === nothing
        @test ctx.lst === l4

        # Modify list
        MinCostFlow.@addstart!(l1, :prv, :nxt, ctx, :fst, :lst)
        MinCostFlow.@addend!(l5, :prv, :nxt, ctx, :fst, :lst)
        MinCostFlow.@remove!(l3, :prv, :nxt, ctx, :fst, :lst)

        @test ctx.fst === l1
        @test l1.prv === nothing
        @test l1.nxt === l2
        @test l2.prv === l1
        @test l2.nxt === l4
        @test l4.prv === l2
        @test l4.nxt === l5
        @test l5.prv === l4
        @test l5.nxt === nothing
        @test ctx.lst === l5

        MinCostFlow.@remove!(l1, :prv, :nxt, ctx, :fst, :lst)

        @test ctx.fst === l2
        @test l2.prv === nothing
        @test l2.nxt === l4
        @test l4.prv === l2
        @test l4.nxt === l5
        @test l5.prv === l4
        @test l5.nxt === nothing
        @test ctx.lst === l5

        MinCostFlow.@remove!(l5, :prv, :nxt, ctx, :fst, :lst)

        @test ctx.fst === l2
        @test l2.prv === nothing
        @test l2.nxt === l4
        @test l4.prv === l2
        @test l4.nxt === nothing
        @test ctx.lst === l4

        MinCostFlow.@remove!(l2, :prv, :nxt, ctx, :fst, :lst)

        @test ctx.fst === l4
        @test l4.prv === nothing
        @test l4.nxt === nothing
        @test ctx.lst === l4

        MinCostFlow.@remove!(l4, :prv, :nxt, ctx, :fst, :lst)

        @test ctx.fst === nothing
        @test ctx.lst === nothing

    end

end
