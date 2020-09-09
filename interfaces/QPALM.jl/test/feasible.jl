using Test

@testset "Feasible" begin

    using QPALM
    using Random
    using LinearAlgebra
    using SparseArrays

    Random.seed!(0)

    n, m = 10, 20

    @testset "Inequalities ($(act) active)" for act in 0:5

        # L(x, y) = (1/2)*(x'*Q*x) + q'*x + y'*(Ax - b)
        # y >= 0 i.e. -min(y, 0) <= tol
        # Ax - b <= 0 i.e. max(Ax - b, 0) <= tol
        # y'*(Ax - b) = 0
        # Q*x + q + A'y = 0

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

        @test results.info.status == :Solved

        @test maximum(-min.(results.y, 0.0)) <= 1e-4
        @test maximum(max.(A * results.x - b, 0.0)) <= 1e-4
        @test abs(dot(results.y, A * results.x - b)) <= 1e-4
        @test norm(Q * results.x + q + A' * results.y, Inf) <= 1e-4

    end

    @testset "Equalities" for _ in 1:5

        # L(x, y) = (1/2)*(x'*Q*x) + q'*x + y'*(Ax - b)
        # Ax - b = 0
        # Q*x + q + A'y = 0

        F = randn(n, n-1)
        Q = sparse(F*F' + 1e-3*I)
        A = sparse(randn(m, n))
        x_star = randn(n)
        y_star = randn(m)
        q = -Q*x_star - A'*y_star
        b = A*x_star

        model = QPALM.Model()
        QPALM.setup!(model, Q=Q, q=q, A=A, bmin=b, bmax=b; Dict{Symbol,Any}(:eps_rel=>0)...)
        results = QPALM.solve!(model)

        @test results.info.status == :Solved
        @test norm(A*results.x - b, Inf) <= 1e-4
        @test norm(Q*results.x + q + A'*results.y, Inf) <= 1e-4

    end

end
