using Test

@testset "Linear programs" begin

    using QPALM
    using Random
    using LinearAlgebra
    using SparseArrays

    Random.seed!(0)

    # Primal LP
    #
    #   minimize    c'x
    #   subject to  Ax = b
    #               x >= 0
    #
    # Dual LP
    #
    #   maximize    b'y
    #   subject to  A'y <= c
    #
    # Optimality conditions
    #
    #   x >= 0              [primal feasibility 1]
    #   Ax = b              [primal feasibility 2]
    #   A'y <= c            [dual feasibility]
    #   X'(c - A'y) = 0     [complementarity slackness]

    function assert_lp_solution(c, A, b, x, y, tol)
        # Check and print solution quality measures (for some reason the
        # returned dual iterate is the negative of the dual LP variable y)

        nonneg = -minimum(min.(0.0, x))
        @test nonneg <= tol

        primal_feasibility = norm(A*x - b)
        @test primal_feasibility <= tol

        dual_feasibility = maximum(max.(0.0, -A'*y - c))
        @test dual_feasibility <= tol

        complementarity = abs(dot(c + A'*y, x))
        @test complementarity <= tol
    end

    T = Float64

    n = 100 # primal dimension
    m = 80 # dual dimension (i.e. number of linear equalities)
    k = 50 # number of active dual constraints (must be 0 <= k <= n)
    x_star = vcat(rand(T, k), zeros(T, n-k)) # primal optimal point
    s_star = vcat(zeros(T, k), rand(T, n-k)) # dual optimal slack variable
    y_star = randn(T, m) # dual optimal point

    A = randn(T, m, n)
    b = A*x_star
    c = A'*y_star + s_star

    @testset "Primal" begin

        A_qpalm = sparse([A; Matrix{T}(I, n, n)])
        bmin = [b; zeros(T, n)]
        bmax = [b; Inf*ones(T, n)]

        model = QPALM.Model()

        QPALM.setup!(model, q=c, A=A_qpalm, bmin=bmin, bmax=bmax, eps_abs=1e-8, eps_rel=1e-8)
        results = QPALM.solve!(model)

        @test results.info.status == :Solved

        assert_lp_solution(c, A, b, results.x, results.y[1:m], 1e-6)

    end

end
