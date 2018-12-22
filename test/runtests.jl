using MinCostFlow
using Test
using LinearAlgebra
using SparseArrays
using MathProgBase
using Clp
using Profile
using Random

include("utils.jl")
include("listutils.jl")

@testset "MinCostFlow" begin

    include("lists.jl")

    @testset "Example networks" begin

        @testset "Bertsekas page 220" begin

            fp = FlowProblem([1,2], [2,3], [5,5], [0,1], [1,0,-1])
            @test MinCostFlow.complementarityslackness(fp) # Initialization should satisfy CS

            solveflows!(fp)
            @test MinCostFlow.complementarityslackness(fp) # Solving should preserve CS
            @test flows(fp) == [1,1]
            @test prices(fp) == [1,1,0]

            lp = linprog(fp)
            @test flows(fp) == lp.flows

        end

        @testset "Bertsekas page 237" begin

            fp = FlowProblem([1,1,2,3,2,3], [2,3,3,2,4,4], [2,2,3,2,1,5],
                              [5,1,4,3,2,0], [3,2,-1,-4])
            @test MinCostFlow.complementarityslackness(fp) # Initialization should satisfy CS

            solveflows!(fp)
            @test MinCostFlow.complementarityslackness(fp) # Solving should preserve CS
            @test flows(fp) == [1,2,2,0,1,3]
            @test prices(fp) == [9,4,0,0]

            lp = linprog(fp)
            @test flows(fp) == lp.flows

        end

        @testset "Ahuja, Magnanti, and Orlin page 421" begin

            fp = FlowProblem([1,1,2,2,2,3,5,4,5], [2,3,3,4,5,5,4,6,6],
                             [8,3,3,7,2,3,4,5,6], [3,2,2,5,2,4,5,3,4],
                             [9,0,0,0,0,-9])
            @test MinCostFlow.complementarityslackness(fp) # Initialization should satisfy CS

            solveflows!(fp)
            @test MinCostFlow.complementarityslackness(fp) # Solving should preserve CS
            @test flows(fp) == [6,3,0,4,2,3,0,4,5]

            lp = linprog(fp)
            @test flows(fp) == lp.flows

        end

    end

    # Note that degeneracy can potentially prevent MinCostFlow results from
    # matching LP results exactly. Instead, we just check that the minimized
    # costs from each method match and that the MinCostFlow solution is in
    # fact feasible.
    @testset "Random Networks" begin
    
        @testset "Previously problematic problems" begin

            fp = FlowProblem([4,3,4,1,1,2,3,4,5,5,5,5], [1,1,2,4,5,5,5,5,1,2,3,4],
                             [19,9,3,9,9999,9999,9999,9999,9999,9999,9999,9999],
                             [1,1,5,2,0,0,0,0,9999,9999,9999,9999], [-2,-7,4,-12,17])
            @test MinCostFlow.complementarityslackness(fp) # Initialization should satisfy CS

            solveflows!(fp)
            @test MinCostFlow.complementarityslackness(fp) # Solving should preserve CS

        end

        N, E = 50, 150

        for i in 1:25

            fp = randomproblem(N, E)
            @test MinCostFlow.complementarityslackness(fp) # Initialization should satisfy CS

            solveflows!(fp)
            @test MinCostFlow.complementarityslackness(fp) # Solving should preserve CS

            lp = linprog(fp)
            @test dot(flows(fp), costs(fp)) == dot(lp.flows, costs(fp))
            @test buildAmatrix(fp) * flows(fp) == .-injections(fp)

        end

    end

    @testset "Random Hotstarts" begin

        N, E = 20, 40

        fp = randomproblem(N, E)
        @test MinCostFlow.complementarityslackness(fp) # Initialization should satisfy CS

        solveflows!(fp)
        @test MinCostFlow.complementarityslackness(fp) # Solving should preserve CS

        lp = linprog(fp)
        @test dot(flows(fp), costs(fp)) == dot(lp.flows, costs(fp))
        @test buildAmatrix(fp) * flows(fp) == .-injections(fp)

        # Randomly modify problem and re-solve
        for i in 1:25

            # Update injections and rebalance at fallback node
            for n in 1:N
                updateinjection!(fp.nodes[n], fp.nodes[N+1],
                                 fp.nodes[n].injection + rand(-3:3))
            end

            # Update flow limits
            for e in 1:E
                updateflowlimit!(fp.edges[e], max(0, fp.edges[e].limit + rand(-3:3)))
            end

            # Ensure updates preserved CS
            @test MinCostFlow.complementarityslackness(fp)

            # Re-solve
            solveflows!(fp)
            @test MinCostFlow.complementarityslackness(fp)

            lp = linprog(fp)
            @test dot(flows(fp), costs(fp)) == dot(lp.flows, costs(fp))
            @test buildAmatrix(fp) * flows(fp) == .-injections(fp)

        end

    end

    if true

        Random.seed!(1234)
        @profile zeros(1)
        Profile.clear()
        N = 40; E = 100
        println("n = $N, e = $E")
        fp = randomproblem(N, E)
        @profile solveflows!(fp)

        for _ in 1:9

            # Update injections and rebalance at fallback node
            for n in 1:N
                updateinjection!(fp.nodes[n], fp.nodes[N+1],
                                 fp.nodes[n].injection + rand(-3:3))
            end

            # Update flow limits
            for e in 1:E
                updateflowlimit!(fp.edges[e], max(0, fp.edges[e].limit + rand(-3:3)))
            end

            @profile solveflows!(fp)

        end

        Profile.print(maxdepth=14)

    end

end

