using MinCostFlow
using Test
using LinearAlgebra
using SparseArrays
using MathProgBase
using Clp

include("utils.jl")

@testset "MinCostFlow" begin

    @testset "Example networks" begin

        @testset "Bertsekas page 220" begin

            fp = FlowProblem([1,2], [2,3], [5,5], [0,1], [1,0,-1])

            solveflows!(fp)
            @test fp.flows == [1,1]
            @test fp.shadowprices == [1,1,0]

            lp = linprog(fp)
            @test fp.flows == lp.flows

        end

        @testset "Bertsekas page 237" begin

            fp = FlowProblem([1,1,2,3,2,3], [2,3,3,2,4,4], [2,2,3,2,1,5],
                              [5,1,4,3,2,0], [3,2,-1,-4])

            solveflows!(fp)
            @test fp.flows == [1,2,2,0,1,3]
            @test fp.shadowprices == [9,4,0,0]

            lp = linprog(fp)
            @test fp.flows == lp.flows

        end

        @testset "Ahuja, Magnanti, and Orlin page 421" begin

            fp = FlowProblem([1,1,2,2,2,3,5,4,5], [2,3,3,4,5,5,4,6,6],
                             [8,3,3,7,2,3,4,5,6], [3,2,2,5,2,4,5,3,4],
                             [9,0,0,0,0,-9])

            solveflows!(fp)
            @test fp.flows == [6,3,0,4,2,3,0,4,5]

            lp = linprog(fp)
            @test fp.flows == lp.flows

        end

    end

    # Note that degeneracy can potentially prevent MinCostFlow results from
    # matching LP results exactly. Instead, we just check that the minimized
    # costs from each method match and that the MinCostFlow solution is in
    # fact feasible.
    @testset "Random Networks" begin
    
        N, E = 50, 150
        fp = randomproblem(N, E)

        solveflows!(fp)
        lp = linprog(fp)
        @test dot(fp.flows, fp.costs) == dot(lp.flows, fp.costs)
        @test buildAmatrix(fp) * fp.flows == .-fp.injections

        # Randomly modify problem and re-solve
        for _ in 1:20

            # Update injections and rebalance at fallback node
            for n in 1:N
                updateinjection!(fp, fp.injections[n] + rand(-3:3), n, N+1)
            end

            # Update flow limits
            for e in 1:E
                updateflowlimit!(fp, max(0, fp.limits[e] + rand(-3:3)), e)
            end

            # Resolve
            solveflows!(fp)
            lp = linprog(fp)
            @test dot(fp.flows, fp.costs) == dot(lp.flows, fp.costs)
            @test buildAmatrix(fp) * fp.flows == .-fp.injections

        end

    end

end
