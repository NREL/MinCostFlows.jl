using Test
using MinCostFlow

@testset "MinCostFlow" begin

    @testset "Bertsekas page 220" begin
        fp = FlowProblem([1,2], [2,3], [5,5], [0,1], [1,0,-1])
        solveflows!(fp)
        @test fp.flows == [1,1]
        @test fp.shadowprices == [1,1,0]
    end

    @testset "Bertsekas page 237" begin
        fp = FlowProblem([1,1,2,3,2,3], [2,3,3,2,4,4], [2,2,3,2,1,5],
                          [5,1,4,3,2,0], [3,2,-1,-4])
        solveflows!(fp)
        @test fp.flows == [1,2,2,0,1,3]
        @test fp.shadowprices == [9,4,0,0]
    end

    @testset "Ahuja, Magnanti, and Orlin page 421" begin
        fp = FlowProblem([1,1,2,2,2,3,5,4,5], [2,3,3,4,5,5,4,6,6],
                         [8,3,3,7,2,3,4,5,6], [3,2,2,5,2,4,5,3,4],
                         [9,0,0,0,0,-9])
        solveflows!(fp)
        @test fp.flows == [6,3,0,4,2,3,0,4,5]
        @test_broken fp.shadowprices == [0,-3,-2,-8,-7,-11]
    end

end
