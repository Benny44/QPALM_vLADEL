using Test

@testset "Update" begin

    using QPALM
    using LinearAlgebra
    using SparseArrays

    q = [0.1, 1.0]
    A = sparse([1.0 0.0; 0.0 1.0; 1.0 1.0])
    bmin = [0.0, 0.0, -Inf]
    bmax = [Inf, Inf, 1.0]

    model = QPALM.Model()

    
    QPALM.setup!(model, q=q, A=A, bmin=bmin, bmax=bmax; Dict{Symbol,Any}(:eps_rel=>0,:eps_abs=>1e-6)...)

    results = QPALM.solve!(model)

    @test isapprox(results.x, [0.0, 0.0], atol=1e-6)

    # @testset "Settings" begin
    #
    #     QPALM.update!(model; max_iter=1042)
    #
    #     settings = unsafe_load(workspace.settings)
    #
    #     @test settings.max_iter == 1042
    #
    # end

    @testset "Bounds" begin

        QPALM.update!(model; bmin=[0.2, 0.0, -Inf])

        results = QPALM.solve!(model)

        @test isapprox(results.x, [0.2, 0.0], atol=1e-6)

    end

    @testset "Linear term" begin

        QPALM.update!(model; q=[-0.1, 1.0])

        results = QPALM.solve!(model)

        @test isapprox(results.x, [1.0, 0.0], atol=1e-6)

    end

end
