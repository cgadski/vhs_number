using Test
using Veronese
using LinearAlgebra: I

@testset verbose = true begin
    @testset "annihilate_subspace" begin
        V = random_subspace(10, 3)
        Vc = annihilate_subspace(V)

        @test size(Vc) == (10, 7)
        @test V' * V ≈ I
        @test Vc' * Vc ≈ I
        @test Vc' * V ≈ zeros(7, 3) atol = 1e-15
    end

    @testset "complement_subspace" begin
        V = random_subspace(10, 3)
        V, Vc = complement_subspace([V V])

        @test size(V) == (10, 3)
        @test size(Vc) == (10, 7)
        @test V' * V ≈ I
        @test Vc' * Vc ≈ I
        @test Vc' * V ≈ zeros(7, 3) atol = 1e-15
    end

    @testset "coordinate_observation" begin
        V = random_subspace(10, 3)
        x = random_vector(3)
        W = coordinate_observation(V * x, random_coordinates(10, 5))
        S = svd(V' * W).S

        @test W' * W ≈ I
        @test S[1] ≈ 1 atol = 1e-15
        @test S[2] < 0.99
    end

    @testset "veronese, un_veronese" begin
        V = random_subspace(10, 3)
        x = random_vector(3)
        y = veronese(V * x)
        x_prime = un_veronese(y, 10)

        @test y ≈ veronese(V) * veronese(x)
        @test abs(x_prime' * V * x) ≈ 1 atol = 1e-15
    end

    @testset "subspace_arrangement" begin
        n, k, d = 10, 3, 3
        subspaces, V, Vc = subspace_arrangement(n, k, d)
        x = random_vector(d)

        @test size(V, 2) == k * binomial(d + 1, 2)
        @test size(Vc, 2) == binomial(n + 1, 2) - k * binomial(d + 1, 2)
        @test norm(Vc' * veronese(subspaces[:, :, 1])) ≈ 0 atol = 1e-15
    end

    @testset "incidence_relations" begin
        V = random_subspace(10, 3)
        x = random_vector(3)
        W = [random_subspace(10, 6) (V * x)]

        Vc = annihilate_subspace(V)
        W, Wc = complement_subspace(W)
        rels = incidence_relations(W, Wc, V, Vc)
        println("relations")
        # @test size(rels) == (7 * 3, 1)
    end
end
