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
            @test MinCostFlow.complementaryslackness(fp) # Initialization should satisfy CS

            lp = linprog(fp)
            solveflows!(fp)

            @test MinCostFlow.complementaryslackness(fp) # Solving should preserve CS
            @test flows(fp) == [1,1]
            @test prices(fp) == [1,1,0]
            @test flows(fp) == lp.flows

        end

        @testset "Bertsekas page 237" begin

            fp = FlowProblem([1,1,2,3,2,3], [2,3,3,2,4,4], [2,2,3,2,1,5],
                              [5,1,4,3,2,0], [3,2,-1,-4])
            @test MinCostFlow.complementaryslackness(fp) # Initialization should satisfy CS

            lp = linprog(fp)
            solveflows!(fp)

            @test MinCostFlow.complementaryslackness(fp) # Solving should preserve CS
            @test flows(fp) == [1,2,2,0,1,3]
            @test prices(fp) == [9,4,0,0]
            @test flows(fp) == lp.flows

        end

        @testset "Ahuja, Magnanti, and Orlin page 421" begin

            fp = FlowProblem([1,1,2,2,2,3,5,4,5], [2,3,3,4,5,5,4,6,6],
                             [8,3,3,7,2,3,4,5,6], [3,2,2,5,2,4,5,3,4],
                             [9,0,0,0,0,-9])
            @test MinCostFlow.complementaryslackness(fp) # Initialization should satisfy CS

            lp = linprog(fp)
            solveflows!(fp)

            @test MinCostFlow.complementaryslackness(fp) # Solving should preserve CS
            @test flows(fp) == [6,3,0,4,2,3,0,4,5]
            @test flows(fp) == lp.flows

        end

    end

    # Note that degeneracy can potentially prevent MinCostFlow results from
    # matching LP results exactly. Instead, we just check that the minimized
    # costs from each method match and that the MinCostFlow solution is in
    # fact feasible.

    @testset "Previously problematic problems" begin

        fp = FlowProblem(
            [7, 5, 23, 36, 14, 2, 27, 6, 4, 1, 15, 27, 37, 11, 4, 21, 18, 37,
             30, 27, 5, 10, 10, 16, 20, 16, 24, 31, 8, 24, 35, 15, 20, 30, 35,
             29, 5, 18, 7, 23, 38, 22, 30, 36, 11, 3, 40, 36, 14, 39, 36, 24,
             19, 27, 3, 35, 30, 8, 36, 32, 7, 15, 9, 31, 7, 19, 28, 16, 22, 11,
             18, 2, 20, 3, 20, 1, 4, 25, 25, 3, 38, 18, 27, 24, 10, 19, 24, 23,
             12, 5, 16, 7, 4, 36, 38, 32, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11,
             12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27,
             28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 41, 41,
             41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41,
             41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41,
             41, 41, 41, 41, 41],
            [27, 1, 14, 6, 24, 20, 11, 18, 1, 8, 30, 5, 10, 13, 17, 15, 20, 36,
             7, 26, 29, 22, 36, 6, 31, 22, 17, 22, 37, 21, 2, 22, 24, 5, 30, 4,
             17, 6, 33, 27, 25, 14, 2, 26, 16, 30, 15, 35, 8, 34, 13, 36, 23,
             40, 14, 9, 15, 31, 37, 8, 1, 19, 38, 3, 34, 39, 2, 13, 38, 3, 16,
             6, 18, 32, 7, 38, 5, 10, 5, 31, 30, 7, 12, 3, 16, 5, 15, 4, 9, 36,
             33, 10, 13, 28, 39, 20, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41,
             41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41,
             41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 1, 2, 3,
             4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21,
             22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37,
             38, 39, 40],
            [8, 4, 13, 9, 1, 18, 9, 8, 16, 12, 3, 2, 18, 12, 18, 8, 11, 14, 4,
             16, 11, 19, 14, 20, 3, 19, 5, 14, 7, 11, 18, 9, 6, 2, 19, 5, 14,
             2, 3, 6, 15, 2, 11, 16, 15, 11, 7, 14, 6, 9, 2, 3, 5, 6, 17, 10,
             8, 8, 4, 7, 13, 10, 18, 17, 12, 10, 7, 14, 17, 7, 10, 16, 17, 3,
             18, 19, 4, 2, 13, 19, 5, 8, 5, 17, 15, 16, 12, 1, 1, 19, 12, 10,
             5, 20, 5, 3, 9999, 9999, 9999, 9999, 9999, 9999, 9999, 9999, 9999,
             9999, 9999, 9999, 9999, 9999, 9999, 9999, 9999, 9999, 9999, 9999,
             9999, 9999, 9999, 9999, 9999, 9999, 9999, 9999, 9999, 9999, 9999,
             9999, 9999, 9999, 9999, 9999, 9999, 9999, 9999, 9999, 9999, 9999,
             9999, 9999, 9999, 9999, 9999, 9999, 9999, 9999, 9999, 9999, 9999,
             9999, 9999, 9999, 9999, 9999, 9999, 9999, 9999, 9999, 9999, 9999,
             9999, 9999, 9999, 9999, 9999, 9999, 9999, 9999, 9999, 9999, 9999,
             9999, 9999, 9999, 9999, 9999],
            [5, 4, 3, 4, 1, 4, 4, 4, 5, 5, 2, 3, 2, 5, 3, 5, 4, 3, 1, 4, 5, 4,
             4, 2, 2, 3, 2, 2, 3, 1, 3, 2, 3, 2, 3, 4, 1, 4, 4, 1, 1, 1, 5, 2,
             2, 4, 4, 2, 4, 2, 1, 3, 2, 3, 4, 2, 3, 3, 1, 5, 4, 1, 4, 1, 3, 3,
             4, 1, 4, 4, 3, 5, 3, 2, 5, 3, 5, 1, 1, 4, 3, 1, 3, 4, 3, 1, 2, 1,
             2, 5, 2, 4, 3, 4, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
             0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
             0, 0, 0, 0, 9999, 9999, 9999, 9999, 9999, 9999, 9999, 9999, 9999,
             9999, 9999, 9999, 9999, 9999, 9999, 9999, 9999, 9999, 9999, 9999,
             9999, 9999, 9999, 9999, 9999, 9999, 9999, 9999, 9999, 9999, 9999,
             9999, 9999, 9999, 9999, 9999, 9999, 9999, 9999, 9999],
            [-3, 5, 13, -10, -19, -17, -11, 3, -13, 1, 3, -18, -6, -15, -19, 20,
             1, -15, 11, -10, 7, -1, -18, 15, 1, -14, -19, 10, -6, -8, 3, 16,
             19, -10, -15, -10, 15, 12, 16, -18, 104])

        @test MinCostFlow.complementaryslackness(fp) # Initialization should satisfy CS

        lp = linprog(fp) # Will complain if infeasible
        solveflows!(fp)

        @test MinCostFlow.complementaryslackness(fp) # Solving should preserve CS
        @test dot(flows(fp), costs(fp)) == dot(lp.flows, costs(fp))
        @test buildAmatrix(fp) * flows(fp) == .-injections(fp)



        fp = FlowProblem([4,3,4,1,1,2,3,4,5,5,5,5], [1,1,2,4,5,5,5,5,1,2,3,4],
                         [19,9,3,9,9999,9999,9999,9999,9999,9999,9999,9999],
                         [1,1,5,2,0,0,0,0,9999,9999,9999,9999], [-2,-7,4,-12,17])
        @test MinCostFlow.complementaryslackness(fp) # Initialization should satisfy CS

        lp = linprog(fp) # Will complain if infeasible
        solveflows!(fp)

        @test MinCostFlow.complementaryslackness(fp) # Solving should preserve CS
        @test dot(flows(fp), costs(fp)) == dot(lp.flows, costs(fp))
        @test buildAmatrix(fp) * flows(fp) == .-injections(fp)

    end

    @testset "Random Networks" begin

        N, E = 40, 100

        rand(1234)
        for i in 1:250

            fp = randomproblem(N, E)
            @test MinCostFlow.complementaryslackness(fp) # Initialization should satisfy CS

            lp = linprog(fp) # Will complain if infeasible
            solveflows!(fp)

            @test MinCostFlow.complementaryslackness(fp) # Solving should preserve CS
            @test dot(flows(fp), costs(fp)) == dot(lp.flows, costs(fp))
            @test buildAmatrix(fp) * flows(fp) == .-injections(fp)

        end

    end

    @testset "Random Hotstarts" begin

        N, E = 40, 100

        fp = randomproblem(N, E)
        @test MinCostFlow.complementaryslackness(fp) # Initialization should satisfy CS

        lp = linprog(fp)
        solveflows!(fp)

        @test MinCostFlow.complementaryslackness(fp) # Solving should preserve CS
        @test dot(flows(fp), costs(fp)) == dot(lp.flows, costs(fp))
        @test buildAmatrix(fp) * flows(fp) == .-injections(fp)

        # Randomly modify problem and re-solve
        for i in 1:250

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
            @test MinCostFlow.complementaryslackness(fp)

            # Re-solve
            lp = linprog(fp)
            solveflows!(fp)

            @test MinCostFlow.complementaryslackness(fp)
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

