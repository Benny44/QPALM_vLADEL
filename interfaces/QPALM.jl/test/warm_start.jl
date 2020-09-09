using Test

@testset "Warm start" begin

    using QPALM
    using LinearAlgebra
    using SparseArrays
    using Random

    Random.seed!(0)

    n, m = 10, 20
    act = 5
    F = randn(n, n-1)
    Q = sparse(F*F' + 1e-3*I)
    A = sparse(randn(m, n))
    x_star = randn(n)
    y_star = [rand(act); zeros(m-act)]
    q = -Q*x_star - A'*y_star
    b = [A[1:act, :]*x_star; A[act+1:end, :]*x_star + rand(m-act)]

    model = QPALM.Model()
    QPALM.setup!(model, Q=Q, q=q, A=A, bmax=b)
    results = QPALM.solve!(model)
    @test results.info.iter > 5

    @testset "Not setup" begin
        model = QPALM.Model()

        @test_throws ErrorException QPALM.warm_start!(model, x_warm_start=results.x)
    end

    @testset "Normal usage (1)" begin
        model = QPALM.Model()

        QPALM.setup!(model, Q=Q, q=q, A=A, bmax=b)
        QPALM.warm_start!(model, x_warm_start=results.x, y_warm_start=results.y)
        res = QPALM.solve!(model)
        @test res.info.status == :Solved
        @test res.info.iter < results.info.iter
    end

    @testset "Normal usage (2)" begin
        model = QPALM.Model()

        QPALM.setup!(model, Q=Q, q=q, A=A, bmax=b)
        QPALM.warm_start!(model, y_warm_start=results.y)
        QPALM.warm_start!(model, x_warm_start=results.x)

        res = QPALM.solve!(model)
        @test res.info.status == :Solved
        @test res.info.iter < results.info.iter
    end

end
