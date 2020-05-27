using InteractiveUtils
using Debugger

@testset "general tests                          " begin

    u = ChebCoeffs(3)
    v = ChebCoeffs([9, 8, 7, 6])

    @test length(u) == 4
    
    @test v[0] == 9
    @test v[1] == 8
    @test v[2] == 7
    @test v[3] == 6

    u[0] = 1; @test u[0] == 1
    u[1] = 2; @test u[1] == 2
    u[2] = 3; @test u[2] == 3
    u[3] = 4; @test u[3] == 4

    a = copy(u)
    @test a[0] == 1
    @test a[1] == 2
    @test a[2] == 3
    @test a[3] == 4

    # find Chebychev series expansion of a function and check its derivative
    P = 21
    y = cos.(π*(0:P)/P)
    f(y) = exp.(y)
    f̂ = ChebCoeffs(FFTW.r2r(f.(y), FFTW.REDFT00)/P)
    f̂[0] /= 2

    @test HelmoltzSolvers._ddy(f̂, Val(:left)) ≈ exp(-1)
    @test HelmoltzSolvers._ddy(f̂, Val(:right)) ≈ exp( 1)

    # u .= 2.0 .* v
    # @test u[0] == 18
    # @test u[1] == 16
    # @test u[2] == 14
    # @test u[3] == 12

    # u .= v
    # @test u[0] == 9
    # @test u[1] == 8
    # @test u[2] == 7
    # @test u[3] == 6

    # B = [5 6 7 8; 9 0 1 2]
    # w = ChebCoeffs(view(B, 1, :))
    # println(w.data)
    # u .= w
    # println(u.data)
    # @test u[0] == 5
    # @test u[1] == 6
    # @test u[2] == 7
    # @test u[3] == 8

    # w[0] = 1
    # w[1] = 2
    # w[2] = 3
    # w[3] = 4
    
    # @test B[1, 1] == 1
    # @test B[1, 2] == 2
    # @test B[1, 3] == 3
    # @test B[1, 4] == 4
end